# FEM MATLAB 脚本第一批代码审查与讲解

## 0. 本批已阅读文件

本批文件包含 20 个 `.m` 脚本：

1. `main_mixed_mech.m`
2. `main_mixed_thermal.m`
3. `buildSystem.m`
4. `assembleSystem.m`
5. `buildBCData.m`
6. `getDofInfo.m`
7. `getElementData.m`
8. `getElementDoFs.m`
9. `findEdgeOwner.m`
10. `assembleNaturalBC_Mech.m`
11. `assembleNaturalBC_Thermal.m`
12. `Edge2_MechTraction.m`
13. `Edge2_Thermal.m`
14. `extractNodalFields.m`
15. `averageFieldToNodes.m`
16. `post_Tri3_Mech.m`
17. `post_Quad4_Mech.m`
18. `post_Quad4_Thermal.m`
19. `getElasticMatrix2D.m`
20. `getElementFaces2D.m`

本报告没有覆盖未上传的函数内部实现，例如 `Tri3_Mech.m`、`Quad4_Mech.m`、`Tri3_Thermal.m`、`Quad4_Thermal.m`、`solveSystem.m`、`postprocessModel.m`、绘图函数等。它们属于后续批次或当前目录外依赖。

## 1. 总体结论

这批代码的总体结构清晰，已经形成了一套统一的 2D FEM 教学框架。主线为：

```matlab
main_mixed_xxx
    -> buildSystem
        -> assembleSystem
            -> getElementData
            -> getElementDoFs
            -> Tri3_xxx / Quad4_xxx
        -> assembleNaturalBC_xxx
            -> findEdgeOwner
            -> Edge2_xxx
    -> buildBCData
    -> solveSystem
    -> extractNodalFields
    -> postprocessModel
    -> averageFieldToNodes
    -> plot...
```

核心思想是用 `model` 结构统一保存几何、网格、材料、物理场、边界条件和载荷。纯力学模型使用：

```matlab
model.fieldNames = {'ux','uy'};
```

纯热传导模型使用：

```matlab
model.fieldNames = {'T'};
```

这使得同一套组装函数可以服务不同物理问题。

## 2. 必须改动或必须注意的问题

### 2.1 `main_mixed_thermal.m` 存在本质边界条件冲突，必须修正

`main_mixed_thermal.m` 中四条边界的定温设置如下：

```matlab
% Left boundary x = 0
T = 0.0

% Right boundary x = 1
T = 0.1

% Bottom boundary y = 0
T = 0.1

% Top boundary y = H
T = 0.1
```

问题出现在左下角和左上角：

- 左下角节点同时属于 left boundary 和 bottom boundary
- 左上角节点同时属于 left boundary 和 top boundary

因此它们会被同时施加两个不同温度：

```matlab
T = 0.0
T = 0.1
```

`buildBCData.m` 中已经有冲突检查逻辑：

```matlab
if isPrescribed(dof)
    if abs(bcData.prescribedValues(dof) - val) > 1e-12
        error('Conflicting essential BC detected at global DoF %d.', dof);
    end
end
```

所以原始 `main_mixed_thermal.m` 在构造 `bcData` 时会报错。

我已经单独生成一个不覆盖原文件的修正版：

```text
main_mixed_thermal_fixed.m
```

修正策略：保留左边界角点为 `T = 0.0`，底边和顶边从第二个节点开始施加 `T = 0.1`。

关键改法：

```matlab
% Bottom boundary y = 0
for i = 1:nx
    nid = node(i,0);
    essBC(end+1,1:2) = {nid, 0.1};
end

% Top boundary y = H
for i = 1:nx
    nid = node(i,ny);
    essBC(end+1,1:2) = {nid, 0.1};
end
```

如果你希望底边和顶边优先，则也可以保留底边、顶边角点，把 left boundary 改成 `j = 1:ny-1`。二者都可以，关键是同一个自由度不能有两个不同定值。

### 2.2 当前批次文件不能单独完整运行，因为依赖函数未上传

主脚本调用了以下函数，但本批文件中没有这些函数：

```matlab
Tri3_Mech.m
Quad4_Mech.m
Tri3_Thermal.m
Quad4_Thermal.m
solveSystem.m
postprocessModel.m
plotMesh.m
plotNodalScalar.m
plotDeformedMeshContour.m
```

这不一定是代码错误。用户已说明这是部分脚本，所以这些依赖很可能在后续批次中。当前结论是：

- 组装层、边界层、后处理辅助层可以单独审查
- 两个 main 文件不能仅靠本批文件直接跑通
- 下一批应优先上传单元文件和求解器文件

### 2.3 `post_Tri3_Thermal.m` 当前缺失

`main_mixed_thermal.m` 使用 Tri3 和 Quad4 混合网格，且调用：

```matlab
results = postprocessModel(model, result.u);
Tele = averageFieldToNodes(model, results, 'T');
qN = averageFieldToNodes(model, results, 'heatFlux');
```

本批只有：

```matlab
post_Quad4_Thermal.m
```

如果 `postprocessModel.m` 需要对 Tri3 thermal 元素调用 `post_Tri3_Thermal.m`，那么热传导后处理会失败。若后续批次已经包含该文件，则没有问题。

### 2.4 `solveSystem` 的接口需要在下一批核对

两个 main 文件都写成：

```matlab
[reaction, solution] = solveSystem(globalMatrix, globalVector, bcData);
```

这表示当前主脚本假设 `solveSystem` 的输入顺序是：

```matlab
K, F, bcData
```

输出顺序是：

```matlab
reaction, solution
```

如果后续上传的 `solveSystem.m` 使用别的接口，例如：

```matlab
[solution, reaction] = solveSystem(K, bcData, F)
```

则主脚本需要同步调整。这个接口必须统一，否则结果变量会被交换。

### 2.5 全局 DoF 映射假设节点编号连续

`getElementDoFs.m` 使用：

```matlab
eleDoFs((a-1)*nNodeDoF+i) = (eleNodeIDs(a)-1)*nNodeDoF + i;
```

这要求节点编号满足：

```matlab
1, 2, 3, ..., numNodes
```

当前两个 main 文件正好满足这个要求。如果以后导入外部网格，节点 ID 不连续，代码需要增加 `nodeID -> row index` 映射表。

## 3. FEM 逻辑顺序讲解

### 3.1 主控层：定义物理问题

相关文件：

```matlab
main_mixed_mech.m
main_mixed_thermal.m
```

主脚本的职责是建立 `model`，然后调用通用 FEM 流程。

纯力学主脚本的物理问题是悬臂梁：

- 长度 `L = 10.0`
- 高度 `H = 1.0`
- 厚度 `t = 1.0`
- 杨氏模量 `E = 210e9`
- 泊松比 `nu = 0.30`
- 右端总向下载荷 `P = 1000.0`
- 左端固支
- 右端施加均布边界 traction

它的 FE 变量是：

```matlab
model.fieldNames = {'ux','uy'};
```

对应每个节点两个自由度：

```matlab
ux, uy
```

纯热传导主脚本的物理问题是矩形区域稳态导热：

- `L = 1.0`
- `H = 2.0`
- `k = 0.6`
- `thickness = 1.0`
- 左、右、上、下边界定温
- 体热源为 0

它的 FE 变量是：

```matlab
model.fieldNames = {'T'};
```

对应每个节点一个自由度：

```matlab
T
```

考试重点：主脚本通常不考整段手写，但很可能考如何定义 `model.fieldNames`、如何生成单元连接表、如何设置本质边界和自然边界。

### 3.2 自由度与单元注册层

相关文件：

```matlab
getDofInfo.m
getElementData.m
getElementDoFs.m
```

#### `getDofInfo.m`

作用：从 `model.fieldNames` 得到每个节点的自由度数量和字段名映射。

核心代码：

```matlab
labels = model.fieldNames(:)';
dofInfo.labels = labels;
dofInfo.nDoFPerNode = numel(labels);

for i = 1:numel(labels)
    dofInfo.map.(labels{i}) = i;
end
```

概念：自由度编号规则。

例如：

```matlab
model.fieldNames = {'ux','uy'};
```

则：

```matlab
ux -> 1
uy -> 2
```

考试重点：可能要求说明每节点自由度数量如何决定整体矩阵大小。

#### `getElementData.m`

作用：根据元素名返回元素节点数、元素函数句柄、所需物理场。

核心代码：

```matlab
case 'Tri3_Mech'
    eleData.nEleNodes = 3;
    eleData.func = @Tri3_Mech;
    eleData.requiredFields = {'ux','uy'};
```

概念：元素注册表。

`assembleSystem` 不直接知道 Tri3 或 Quad4 的公式，它只通过 `eleData.func` 调用相应单元函数。

兼容性检查：

```matlab
if ~isequal(eleData.requiredFields, dofInfo.labels)
    error('Element %s is incompatible with model.fieldNames.', eleName);
end
```

这能防止把热单元放进纯力学模型，或把力学单元放进热模型。

考试重点：很可能考元素类型、节点数、自由度数、函数句柄之间的关系。

#### `getElementDoFs.m`

作用：把一个单元的节点编号转换为整体自由度编号。

核心代码：

```matlab
for a = 1:nEleNodes
    for i = 1:nNodeDoF
        eleDoFs((a-1)*nNodeDoF+i) = (eleNodeIDs(a)-1)*nNodeDoF + i;
    end
end
```

例子：二维力学单元，节点 `[3, 7, 8]`，每节点 2 个自由度，则全局 DoF 为：

```matlab
[5, 6, 13, 14, 15, 16]
```

物理意义：局部矩阵的行列要装配到整体矩阵的对应位置。

考试重点：这是高频代码段。考试可能让你补全自由度映射循环。

### 3.3 组装层：从单元矩阵到整体矩阵

相关文件：

```matlab
buildSystem.m
assembleSystem.m
```

#### `assembleSystem.m`

作用：遍历所有单元集，调用单元函数，装配全局矩阵和向量。

核心 FEM 公式：

```matlab
K(eleDoFs, eleDoFs) = K(eleDoFs, eleDoFs) + Ke;
F(eleDoFs) = F(eleDoFs) + Fe;
```

代码中对应：

```matlab
globalMatrix(eleDoFs, eleDoFs) = globalMatrix(eleDoFs, eleDoFs) + eleMatrix;
globalVector(eleDoFs) = globalVector(eleDoFs) + eleVector;
```

这就是有限元 assembly 的核心。

完整逻辑：

1. 读取节点和单元集
2. 检查所有单元集自由度数相同
3. 初始化整体矩阵 `globalMatrix`
4. 初始化整体向量 `globalVector`
5. 遍历每个元素集
6. 读取当前单元的节点编号
7. 根据节点编号取坐标
8. 根据 `model.section` 取材料
9. 根据 `model.eleLoadData` 取体载荷或热源
10. 调用单元函数
11. 得到 `eleDoFs`
12. 装配到整体系统
13. 最后叠加 `model.nodalVector`

考试重点：装配代码几乎一定是重点，尤其是 `eleDoFs` 和矩阵加和。

#### `buildSystem.m`

作用：在域内单元组装结果上继续加入自然边界条件。

核心代码：

```matlab
[globalMatrix, globalVector] = assembleSystem(model);

if isequal(fieldNames, {'T'})
    [bcMatrix, bcVector] = assembleNaturalBC_Thermal(model);
    globalMatrix = globalMatrix + bcMatrix;
    globalVector = globalVector + bcVector;

elseif isequal(fieldNames, {'ux','uy'})
    bcVector = assembleNaturalBC_Mech(model);
    globalVector = globalVector + bcVector;
end
```

物理意义：

- 域内积分贡献来自单元内部
- 自然边界贡献来自边界积分
- 本质边界不在这里施加

机械自然边界只加向量，因为 traction 是外力项。

热自然边界可加矩阵和向量，因为 Robin 边界通常写成：

```text
M T = S
```

考试重点：会区分本质边界和自然边界。

### 3.4 本质边界条件处理

相关文件：

```matlab
buildBCData.m
```

作用：把 `model.essBC` 变成求解器可以用的已知自由度信息。

当前支持格式：

```matlab
model.essBC = {
    [node list], [value vector]
};
```

机械例子：

```matlab
model.essBC = {
    leftNodes, [0,0]
};
```

热学例子：

```matlab
essBC(end+1,1:2) = {nid, 0.1};
```

核心代码：

```matlab
dof = (nodeID-1)*nNodeDoF + j;
val = bcVals(j);
```

如果同一自由度重复施加相同值，允许。如果重复施加不同值，报错。

考试重点：

- Dirichlet 条件对应已知主变量
- 力学中是位移约束
- 热传导中是温度约束
- 已知自由度需要从未知方程中分离出来求解

### 3.5 自然边界条件处理

相关文件：

```matlab
assembleNaturalBC_Mech.m
assembleNaturalBC_Thermal.m
findEdgeOwner.m
Edge2_MechTraction.m
Edge2_Thermal.m
```

#### `findEdgeOwner.m`

作用：检查给定两个节点是否真的构成某个单元的边。

对于 Tri3：

```matlab
[1 2]
[2 3]
[3 1]
```

对于 Quad4：

```matlab
[1 2]
[2 3]
[3 4]
[4 1]
```

这一步能防止把不相邻节点误当成边界边。

考试重点：边界边必须来自单元边，不能随便选两个节点。

#### `Edge2_MechTraction.m`

作用：2 节点线边界单元，把均布 traction 转成等效节点力。

对应公式：

```text
fe = ∫ N^T tbar thickness dGamma
```

代码逻辑：

```matlab
N = [(1-s)/2, (1+s)/2];
J = L/2;
Nmat = [N1 0 N2 0;
        0 N1 0 N2];
edgeVector = edgeVector + Nmat' * [tx;ty] * J * w * thickness;
```

对于常值 traction，两个端点分到的力各占一半。

考试重点：可能考 2 点 Gauss 线积分和等效节点力。

#### `assembleNaturalBC_Mech.m`

作用：把边界链拆成多个 2 节点边界单元，并装配到整体力向量。

核心代码：

```matlab
edgeDoFs = [(nodeA-1)*nNodeDoF+1;
            (nodeA-1)*nNodeDoF+2;
            (nodeB-1)*nNodeDoF+1;
            (nodeB-1)*nNodeDoF+2];

bcVector(edgeDoFs) = bcVector(edgeDoFs) + edgeVector;
```

概念：边界积分也要做 assembly。

#### `Edge2_Thermal.m`

作用：2 节点热边界单元，返回边界矩阵和边界向量。

代码中采用：

```matlab
edgeMatrix = edgeMatrix + (N' * M * N) * J * w * thickness;
edgeVector = edgeVector + (N' * S) * J * w * thickness;
```

适合表达 Robin 或广义自然边界。

考试重点：要知道热边界可能影响矩阵，也可能只影响右端项，取决于边界形式。

#### `assembleNaturalBC_Thermal.m`

作用：和机械自然边界类似，但需要同时装配矩阵和向量。

核心代码：

```matlab
bcMatrix(edgeDoFs, edgeDoFs) = bcMatrix(edgeDoFs, edgeDoFs) + edgeMatrix;
bcVector(edgeDoFs) = bcVector(edgeDoFs) + edgeVector;
```

考试重点：自然边界条件的装配方式与域内单元矩阵装配方式相同。

### 3.6 求解层

相关文件：

```matlab
solveSystem.m
```

该文件本批未上传，但从主脚本可推断它需要完成：

1. 读取已知自由度 `bcData.prescribedDoFs`
2. 找到未知自由度
3. 分块矩阵
4. 求解未知位移或温度
5. 回代计算完整反力或热反力

标准分块形式：

```text
[Kuu Kuk] [u_unknown] = [F_unknown]
[Kku Kkk] [u_known  ]   [F_known  ]
```

未知量求解：

```text
u_unknown = Kuu \ (F_unknown - Kuk*u_known)
```

考试重点：这是极高频内容。笔试很可能给出矩阵分块，让你补代码或解释已知位移自由度如何处理。

### 3.7 后处理层

相关文件：

```matlab
extractNodalFields.m
averageFieldToNodes.m
post_Tri3_Mech.m
post_Quad4_Mech.m
post_Quad4_Thermal.m
```

#### `extractNodalFields.m`

作用：把整体解向量拆成按节点排列的场变量。

对于力学：

```matlab
nodal.ux
nodal.uy
```

对于热学：

```matlab
nodal.T
```

核心代码：

```matlab
base = (n-1)*nNodeDoF;
nodal.(fieldNames{i})(n) = solution(base+i);
```

考试重点：会从整体解向量恢复节点物理量。

#### `averageFieldToNodes.m`

作用：把单元中心或单元积分点结果平均到节点。

例如单元应力是元素结果，但绘制云图常常需要节点值。代码做法：

```matlab
sumVals(idx,:) = sumVals(idx,:) + val;
countVals(idx) = countVals(idx) + 1;
```

最终：

```matlab
nodalField(i,:) = sumVals(i,:) / countVals(i);
```

物理意义：应力、热流等派生量通常在单元内计算，再通过平均得到节点可视化结果。

考试重点：后处理平均概念可能会问，但手写整段概率低。

#### `post_Tri3_Mech.m`

作用：计算 Tri3 机械单元的应变、应力和 von Mises 应力。

Tri3 位移插值为线性，单元内应变为常数。

核心 B 矩阵：

```matlab
B = [b1 0 b2 0 b3 0;
     0 c1 0 c2 0 c3;
     c1 b1 c2 b2 c3 b3] / detA;
```

对应公式：

```text
strain = B * eleSol
stress = D * strain
```

von Mises：

```matlab
sqrt(stress(1)^2 - stress(1)*stress(2) + stress(2)^2 + 3*stress(3)^2)
```

考试重点：Tri3 的 B 矩阵、常应变性质、`stress = D*strain` 是高频。

#### `post_Quad4_Mech.m`

作用：在 Quad4 单元中心计算应变、应力和 von Mises。

Quad4 形函数在自然坐标中定义，先求：

```matlab
dNds, dNdt
```

然后建立 Jacobian：

```matlab
J = [dNds'*x, dNds'*y;
     dNdt'*x, dNdt'*y];
```

再把自然坐标导数转成实际坐标导数：

```matlab
gradN = inv(J) * [dNds'; dNdt'];
```

最后构造 B 矩阵：

```matlab
B = [dNdx1 0 ...;
     0 dNdy1 ...;
     dNdy1 dNdx1 ...];
```

考试重点：Q4 的 Jacobian、导数变换、B 矩阵构造是高频。

#### `post_Quad4_Thermal.m`

作用：在 Quad4 中心计算温度、温度梯度和热流。

核心公式：

```text
gradT = gradN * eleSol
q = -Kcond * gradT
```

这对应 Fourier 导热定律：

```text
q = -k grad T
```

考试重点：热传导中 `gradT` 和 `heatFlux` 的关系非常重要。

### 3.8 材料矩阵层

相关文件：

```matlab
getElasticMatrix2D.m
```

作用：根据 `planeStress` 或 `planeStrain` 返回二维线弹性矩阵 `D`。

平面应力：

```matlab
D = E/(1 - nu^2) * [1   nu  0;
                    nu  1   0;
                    0   0   (1 - nu)/2];
```

平面应变：

```matlab
D = E/((1 + nu)*(1 - 2*nu)) * [1 - nu   nu       0;
                               nu       1 - nu   0;
                               0        0        (1 - 2*nu)/2];
```

考试重点：必须会区分平面应力和平面应变。薄板通常用平面应力，长厚结构截面常用平面应变。

### 3.9 绘图辅助层

相关文件：

```matlab
getElementFaces2D.m
```

作用：根据单元类型返回绘图用的面节点编号。

Tri3：

```matlab
face = eleNodeIDs(1:3);
```

Quad4：

```matlab
face = eleNodeIDs(1:4);
```

概念：数值计算用的是单元连接表，绘图需要面片连接表。

## 4. 两个 main 文件的工作逻辑

### 4.1 `main_mixed_mech.m`

逻辑顺序：

1. 清空环境
2. 定义梁尺寸、材料参数和载荷
3. 定义网格密度 `nx, ny`
4. 定义 Tri3 和 Quad4 单元类型
5. 生成规则矩形节点
6. 将每个网格 cell 转成 Tri3 或 Quad4 单元
7. 写入 `model.elesInfo`
8. 写入材料数据 `model.section`
9. 写入体力 `model.eleLoadData`
10. 左边界固支
11. 右边界施加均布 traction
12. 调用 `buildSystem`
13. 调用 `buildBCData`
14. 调用 `solveSystem`
15. 提取节点位移
16. 后处理应力和 von Mises
17. 与梁理论位移和反力做对比
18. 绘图

关键 FEM 公式：

```text
K u = F
```

机械单元公式：

```text
Ke = ∫ B^T D B t dA
```

自然边界 traction：

```text
fe = ∫ N^T tbar t dGamma
```

考试重点：

- 悬臂梁边界条件设置
- traction 如何由总载荷转换
- reaction check 的物理意义
- FE 位移与 Euler Bernoulli / Timoshenko 梁理论对比

当前该文件本身未发现明显逻辑错误。它的可运行性依赖未上传的单元、求解、后处理和绘图函数。

### 4.2 `main_mixed_thermal.m`

逻辑顺序：

1. 清空环境
2. 定义矩形导热区域
3. 定义导热系数 `k`
4. 生成混合 Tri3 / Quad4 网格
5. 设置每个单元的 `[kxx, kyy]`
6. 设置热源 `heatSource = 0.0`
7. 设置温度本质边界
8. 调用 `buildSystem`
9. 调用 `buildBCData`
10. 调用 `solveSystem`
11. 提取节点温度
12. 后处理热流
13. 绘制温度和热流图

热传导单元公式：

```text
Ke = ∫ B_T^T k B_T t dA
```

热流：

```text
q = -k gradT
```

当前必须修正的问题是左侧角点温度冲突。已提供 `main_mixed_thermal_fixed.m`。

## 5. 考试最可能涉及的关键代码段

### 5.1 自由度映射

高频程度：极高。

```matlab
eleDoFs((a-1)*nNodeDoF+i) = (eleNodeIDs(a)-1)*nNodeDoF + i;
```

必须理解：节点编号转换为整体矩阵行列编号。

### 5.2 单元装配

高频程度：极高。

```matlab
globalMatrix(eleDoFs, eleDoFs) = globalMatrix(eleDoFs, eleDoFs) + eleMatrix;
globalVector(eleDoFs) = globalVector(eleDoFs) + eleVector;
```

必须理解：共享节点处不同单元的刚度贡献会相加。

### 5.3 本质边界条件

高频程度：高。

```matlab
if isPrescribed(dof)
    if abs(oldVal - newVal) > 1e-12
        error(...)
    end
end
```

必须理解：同一自由度不能同时规定两个不同值。

### 5.4 Tri3 B 矩阵

高频程度：极高。

```matlab
B = [b1 0 b2 0 b3 0;
     0 c1 0 c2 0 c3;
     c1 b1 c2 b2 c3 b3] / detA;
```

必须理解：Tri3 是常应变三角形。

### 5.5 Q4 Jacobian 和导数变换

高频程度：极高。

```matlab
J = [dNds'*x, dNds'*y;
     dNdt'*x, dNdt'*y];

gradN = inv(J) * [dNds'; dNdt'];
```

必须理解：自然坐标导数必须通过 Jacobian 转换成物理坐标导数。

### 5.6 机械材料矩阵

高频程度：高。

```matlab
D = getElasticMatrix2D(E, nu, mechType);
stress = D * strain;
```

必须理解：平面应力和平面应变的区别。

### 5.7 热流计算

高频程度：高。

```matlab
gradT = gradN * eleSol;
q = -Kcond * gradT;
```

必须理解：热流方向与温度梯度方向相反。

### 5.8 自然边界边界单元

高频程度：中高。

```matlab
N = [(1-s)/2, (1+s)/2];
J = L/2;
edgeVector = edgeVector + Nmat' * traction * J * w * thickness;
```

必须理解：边界 traction 或热边界条件也需要积分并装配。

## 6. 建议下一批优先上传的文件

为了完整闭合本批代码，建议下一批优先上传：

1. `Tri3_Mech.m`
2. `Quad4_Mech.m`
3. `Tri3_Thermal.m`
4. `Quad4_Thermal.m`
5. `solveSystem.m`
6. `postprocessModel.m`
7. `post_Tri3_Thermal.m`
8. `plotMesh.m`
9. `plotNodalScalar.m`
10. `plotDeformedMeshContour.m`

优先级最高的是前 6 个，因为它们决定主程序是否能完成求解流程。

## 7. 当前保留不改的内容

以下内容没有明显错误，建议暂时不改：

1. `assembleSystem.m` 的整体装配逻辑
2. `getElementDoFs.m` 的自由度映射，前提是节点编号连续
3. `getElementData.m` 的物理场兼容性检查
4. `buildSystem.m` 的按 `fieldNames` 分流逻辑
5. `assembleNaturalBC_Mech.m` 的边界链拆分和装配逻辑
6. `assembleNaturalBC_Thermal.m` 的热边界矩阵和向量装配逻辑
7. `Edge2_MechTraction.m` 的 2 点 Gauss 线积分
8. `Edge2_Thermal.m` 的 2 点 Gauss 线积分
9. `post_Tri3_Mech.m` 的常应变后处理
10. `post_Quad4_Mech.m` 和 `post_Quad4_Thermal.m` 的中心点后处理

## 8. 最简记忆链条

这批代码对应的 FEM 记忆链条如下：

```text
model 定义问题
-> getDofInfo 定义自由度
-> getElementData 选择元素函数
-> assembleSystem 组装域内矩阵
-> assembleNaturalBC 加自然边界
-> buildBCData 处理本质边界
-> solveSystem 解 Ku = F
-> extractNodalFields 拆解节点结果
-> post_xxx 求应变、应力、热流
-> averageFieldToNodes 转成节点云图
```

力学核心链条：

```text
u -> strain = B u -> stress = D strain -> K = ∫ B^T D B dΩ
```

热学核心链条：

```text
T -> gradT = B_T T -> q = -k gradT -> K = ∫ B_T^T k B_T dΩ
```

