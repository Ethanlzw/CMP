# FEM MATLAB 脚本第三批代码审查与讲解：热力耦合 Thermomechanical

## 0. 本批已阅读文件

本批上传的是热力耦合 FEM 脚本，主要文件如下：

1. `ReadMe.txt`
2. `main_mixed_thermomech.m`
3. `getElementData.m`
4. `getElementDoFs.m`
5. `getElementFaces2D.m`
6. `Quad4_ThermoMech.m`
7. `post_Tri3_ThermoMech.m`
8. `post_Quad4_ThermoMech.m`
9. `postprocessModel.m`
10. `solveSystem.m`
11. `assembleNaturalBC_ThermoMech.m`
12. `assembleSystem.m`
13. `averageFieldToNodes.m`
14. `buildBCData.m`
15. `buildSystem.m`
16. `Edge2_ThermoMech.m`
17. `extractNodalFields.m`
18. `findEdgeOwner.m`
19. `getDofInfo.m`
20. `getElasticMatrix2D.m`

另外，本批代码中 `getElementData.m` 和主程序都引用了：

```matlab
Tri3_ThermoMech.m
```

但当前文件夹没有上传这个单元函数。当前默认 `nTriCols = 0`，所以主程序实际只使用 Quad4 单元，可以暂时避开 Tri3 单元计算。但如果要真正使用 mixed mesh，也就是 Tri3 + Quad4 混合网格，必须补齐 `Tri3_ThermoMech.m`。

---

## 1. 总体结论

这套热力耦合代码的总体框架是正确的。它实现的是**线性稳态一向热力耦合**：

```text
温度场 T 影响机械场 ux, uy
机械变形不反过来影响稳态热传导
```

对应的整体方程是：

```text
[ Kuu  KuT ] [ u ] = [ fu ]
[  0   KTT ] [ T ]   [ fT ]
```

其中：

```text
Kuu : 机械刚度矩阵
KTT : 热传导矩阵
KuT : 温度引起的等效热应变耦合矩阵
```

这与课程中热力耦合 FEM 的矩阵形式一致。

当前代码有 2 个必须关注的问题和 2 个建议修正点：

### 必须关注的问题 1：缺少 `Tri3_ThermoMech.m`

若 `nTriCols > 0`，程序会生成三角形热力耦合单元，但没有对应的 `Tri3_ThermoMech.m`，会导致运行失败。

### 必须关注的问题 2：`Quad4_ThermoMech.m` 没有把 `Tref` 加入单元右端项

代码中后处理已经使用：

```matlab
strainTh = alpha * (Tc - Tref) * [1; 1; 0];
```

说明物理模型采用的是：

```text
thermal strain = alpha * (T - Tref) * [1; 1; 0]
```

但单元矩阵中只写了 `KuT * T`，没有加入 reference temperature 对右端项的贡献。

当 `Tref = 0` 时，原代码不会出数值问题。

当 `Tref ~= 0` 时，原代码的单元方程和后处理物理意义不一致。

修正方法是在 Gauss 积分中加入：

```matlab
fMech = fMech - (Bu' * D * epsTh) * Tref * detJ * w * thickness;
```

### 建议修正点 1：Q4 形函数节点顺序建议与主程序统一

主程序中 Q4 单元节点顺序为：

```matlab
[n1, n2, n3, n4] = [bottom-left, bottom-right, top-right, top-left]
```

课程和多数 FEM 程序的标准自然坐标对应为：

```text
node 1 : (-1,-1)
node 2 : ( 1,-1)
node 3 : ( 1, 1)
node 4 : (-1, 1)
```

原始 `Quad4_ThermoMech.m` 使用的是旋转后的自然坐标定义。在规则矩形单元中通常还能得到可用结果，因为映射相当于自然坐标整体旋转。但为了和课程公式、主程序节点顺序、前面纯力学/纯热传导脚本保持一致，建议统一到标准 Q4 写法。

### 建议修正点 2：`assembleSystem.m` 应先跳过空单元集，再解析单元函数

原始代码中：

```matlab
eleData = getElementData(eleName, model);
...
if isempty(elesSet)
    continue;
end
```

建议改成先判断 `elesSet` 是否为空，再调用 `getElementData`。这样当 `elesSet1` 为空且 `Tri3_ThermoMech.m` 暂时不存在时，不会因为解析空单元集而产生潜在问题。

---

## 2. 热力耦合 FEM 的物理与数学逻辑

## 2.1 主变量

热力耦合二维问题每个节点有 3 个自由度：

```text
ux, uy, T
```

因此：

```matlab
model.fieldNames = {'ux','uy','T'};
```

每个节点的局部自由度顺序是：

```text
[ux, uy, T]
```

一个 Tri3 单元有 3 个节点，总自由度数为：

```text
3 nodes * 3 dofs/node = 9 dofs
```

一个 Quad4 单元有 4 个节点，总自由度数为：

```text
4 nodes * 3 dofs/node = 12 dofs
```

---

## 2.2 热力耦合的本构关系

二维小变形应变为：

```text
strain = [exx, eyy, gamma_xy]^T
```

机械应变来自位移：

```text
strain = Bu * u
```

热应变为：

```text
strainTh = alpha * (T - Tref) * [1; 1; 0]
```

应力为：

```text
stress = D * (strain - strainTh)
```

也就是：

```text
stress = D * (Bu*u - alpha*(T-Tref)*[1;1;0])
```

热传导关系为：

```text
q = -Kcond * gradT
```

其中：

```text
Kcond = [kxx 0; 0 kyy]
```

---

## 2.3 单元矩阵结构

热力耦合单元矩阵应写成：

```text
[ Kuu  KuT ]
[  0   KTT ]
```

其中：

```text
Kuu = integral(Bu' * D * Bu dOmega)
```

```text
KTT = integral(BT' * Kcond * BT dOmega)
```

```text
KuT = - integral(Bu' * D * alpha*[1;1;0] * NT dOmega)
```

这里负号来自：

```text
stress = D * (mechanical strain - thermal strain)
```

这段是考试很可能考的核心公式。

---

## 3. 文件逐个审查与讲解

## 3.1 `ReadMe.txt`

### 作用

说明这套代码是 SDU Computational Multiphysics 课程中用于线性 FEA 教学的 MATLAB 代码，并列出了需要学生完成的模板文件。

### 关键点

ReadMe 中列出的模板包括：

```text
getElasticMatrix2D.m
Edge2_ThermoMech.m
assembleNatureBC_ThermoMech.m
assmbleSystem.m
Quad4_ThermoMech.m
post_Quad4_ThermoMech.m
solveSystem.m
```

其中有两个拼写问题需要注意：

```text
assembleNatureBC_ThermoMech.m
assmbleSystem.m
```

实际代码文件名是：

```text
assembleNaturalBC_ThermoMech.m
assembleSystem.m
```

这属于 ReadMe 拼写问题，不是代码错误。

---

## 3.2 `main_mixed_thermomech.m`

### 作用

这是热力耦合主程序，用于做一个**自由热膨胀 patch test**。

它定义：

```text
几何尺寸 L, H
材料参数 E, nu, alpha, kxx, kyy
参考温度 Tref
温升 DeltaT
混合网格参数 nx, ny, nTriCols
边界精确位移和温度
```

### 工作逻辑

主流程为：

```text
1. 定义 model.fieldNames = {'ux','uy','T'}
2. 定义 Tri3_ThermoMech 和 Quad4_ThermoMech 两类单元
3. 生成规则矩形节点
4. 根据 nTriCols 生成 Tri3 + Quad4 混合网格
5. 设置材料 model.section
6. 设置单元载荷 model.eleLoadData
7. 给所有边界节点施加精确自由热膨胀解
8. buildSystem
9. buildBCData
10. solveSystem
11. postprocessModel
12. 与精确解比较误差
13. 绘图
```

### 为什么这个例子是自由热膨胀验证

精确解为：

```text
T = Tref + DeltaT
ux = alpha * DeltaT * x
uy = alpha * DeltaT * y
```

对应机械应变为：

```text
exx = alpha * DeltaT
eyy = alpha * DeltaT
gamma_xy = 0
```

热应变也是：

```text
thermal strain = [alpha*DeltaT, alpha*DeltaT, 0]^T
```

所以：

```text
stress = D * (strain - thermal strain) = 0
```

因此这个验证的目标是：

```text
位移误差接近 0
温度误差接近 0
应变等于热应变
应力接近 0
von Mises 接近 0
```

### 当前限制

原始代码中：

```matlab
nTriCols = 0;
```

所以默认只使用 Quad4 单元。

如果改成：

```matlab
nTriCols > 0
```

就必须提供 `Tri3_ThermoMech.m`。

### 考试重点

这段最可能考：

```matlab
ux = alpha * DeltaT * x;
uy = alpha * DeltaT * y;
essBC(end+1,1:2) = {nid, [ux, uy, Tuni]};
```

它体现了自由热膨胀的精确解。

---

## 3.3 `getElementData.m`

### 作用

这是单元注册表。它根据单元名称返回：

```text
单元节点数
单元函数句柄
所需自由度字段
单元自由度数
```

### 热力耦合部分

```matlab
case 'Tri3_ThermoMech'
    eleData.nEleNodes = 3;
    eleData.func = @Tri3_ThermoMech;
    eleData.requiredFields = {'ux','uy','T'};

case 'Quad4_ThermoMech'
    eleData.nEleNodes = 4;
    eleData.func = @Quad4_ThermoMech;
    eleData.requiredFields = {'ux','uy','T'};
```

### FEM 意义

这保证了热力耦合单元只能在：

```matlab
model.fieldNames = {'ux','uy','T'};
```

的模型中使用。

### 代码评价

逻辑正确。

但当前文件夹缺少 `Tri3_ThermoMech.m`，因此 mixed mesh 不完整。

### 考试重点

可能考单元注册和自由度兼容检查：

```matlab
if ~isequal(eleData.requiredFields, dofInfo.labels)
    error('Element %s is incompatible with model.fieldNames.', eleName);
end
```

---

## 3.4 `getElementDoFs.m`

### 作用

把单元节点编号转成全局自由度编号。

核心代码：

```matlab
for a = 1:nEleNodes
    for i = 1:nNodeDoF
        eleDoFs((a-1)*nNodeDoF+i) = (eleNodeIDs(a)-1)*nNodeDoF + i;
    end
end
```

### 热力耦合例子

若节点 5 的自由度顺序是：

```text
ux5, uy5, T5
```

则它的全局自由度为：

```text
(5-1)*3+1 = 13 -> ux5
(5-1)*3+2 = 14 -> uy5
(5-1)*3+3 = 15 -> T5
```

### 代码评价

正确。

### 考试重点

这是最可能手写的代码之一。

---

## 3.5 `getElementFaces2D.m`

### 作用

返回 2D 单元的面连接关系，用于绘图。

```matlab
Tri3 -> nodeIDs(1:3)
Quad4 -> nodeIDs(1:4)
```

### 代码评价

正确。

---

## 3.6 `getElasticMatrix2D.m`

### 作用

返回二维线弹性材料矩阵 `D`。

### 平面应力

```matlab
D = E/(1 - nu^2) * [1   nu  0;
                    nu  1   0;
                    0   0   (1 - nu)/2];
```

适合薄板。

### 平面应变

```matlab
D = E/((1 + nu)*(1 - 2*nu)) * [1 - nu   nu       0;
                               nu       1 - nu   0;
                               0        0        (1 - 2*nu)/2];
```

适合长厚结构截面。

### 代码评价

正确。

### 考试重点

这个文件是 ReadMe 明确列出的模板文件之一，也是考试可能要求补全的函数。

---

## 3.7 `Quad4_ThermoMech.m`

### 作用

这是 Q4 四节点热力耦合单元函数。它输出：

```text
eleMatrix : 12 x 12 单元矩阵
eleVector : 12 x 1 单元载荷向量
info      : 单元信息
```

### 自由度排序

每个节点：

```text
[ux, uy, T]
```

四节点单元总自由度顺序：

```text
[ux1 uy1 T1 ux2 uy2 T2 ux3 uy3 T3 ux4 uy4 T4]
```

因此机械自由度为：

```matlab
mechDofs = [1 2 4 5 7 8 10 11];
```

温度自由度为：

```matlab
tempDofs = [3 6 9 12];
```

### 核心矩阵

```matlab
Kuu = Kuu + (Bu' * D * Bu) * detJ * weights(igp) * thickness;
KTT = KTT + (BT' * Kcond * BT) * detJ * weights(igp) * thickness;
KuT = KuT - (Bu' * D * epsTh) * N' * detJ * weights(igp) * thickness;
```

其中：

```text
Bu : 应变-位移矩阵
BT : 温度梯度矩阵
D  : 弹性矩阵
Kcond : 导热矩阵
N  : 温度插值形函数
```

### 关键物理含义

`KuT` 的含义是：

```text
温度变化 -> 热应变 -> 等效机械载荷
```

负号来自：

```text
stress = D * (strain - thermal strain)
```

### 当前问题

原始文件中有：

```matlab
%fMech = fMech - (Bu' * D * epsTh) * Tref * detJ * weights(igp) * thickness;
```

这行被注释了。

若 `Tref = 0`，结果不受影响。

若 `Tref ~= 0`，应取消注释并保留负号。

修正版中已加入：

```matlab
fMech = fMech - (Bu' * D * epsTh) * Tref * detJ * w * thickness;
```

### 考试重点

最可能考这三行：

```matlab
Kuu = Kuu + (Bu' * D * Bu) * detJ * w * thickness;
KTT = KTT + (BT' * Kcond * BT) * detJ * w * thickness;
KuT = KuT - (Bu' * D * epsTh) * N' * detJ * w * thickness;
```

---

## 3.8 `post_Quad4_ThermoMech.m`

### 作用

对 Q4 热力耦合单元在中心点做后处理。

输出包括：

```text
centroid
T
strain
thermalStrain
stress
vonMises
gradT
heatFlux
```

### 核心后处理公式

```matlab
strain = Bu * u;
Tc = N' * T;
strainTh = alpha * (Tc - Tref) * [1; 1; 0];
stress = D * (strain - strainTh);

gradT = BT * T;
q = -Kcond * gradT;
```

### 代码评价

物理公式正确。

建议与 `Quad4_ThermoMech.m` 使用完全相同的 Q4 形函数顺序，避免单元计算和后处理不一致。

---

## 3.9 `post_Tri3_ThermoMech.m`

### 作用

对 Tri3 热力耦合单元做后处理。

Tri3 是线性三角形单元，所以：

```text
Bu 常数
BT 常数
strain 常数
gradT 常数
stress 常数
heatFlux 常数
```

### 核心代码

```matlab
strain = Bu * u;
Tc = mean(T);
strainTh = alpha * (Tc - Tref) * [1;1;0];
stress = D * (strain - strainTh);

gradT = BT * T;
q = -Kcond * gradT;
```

### 代码评价

后处理逻辑正确。

但由于当前缺少 `Tri3_ThermoMech.m`，这个后处理函数暂时没有对应的单元计算函数。

---

## 3.10 `Edge2_ThermoMech.m`

### 作用

二维边界二节点单元，用于装配热力耦合自然边界条件。

输入：

```matlab
nodeCoords : 2 x 2 [x y]
bcVals     : [tx, ty, M, S]
thickness  : thickness
```

输出：

```matlab
edgeMatrix : 6 x 6
edgeVector : 6 x 1
```

边界自由度顺序：

```text
[ux1 uy1 T1 ux2 uy2 T2]
```

### 机械边界力

```matlab
Nmat_u = [N(1) 0    N(2) 0;
          0    N(1) 0    N(2)];
fu = Nmat_u' * [tx; ty] * J * w(igp) * thickness;
```

这表示边界 traction 被转换成等效节点力。

### 热边界项

```matlab
KT = (N' * M * N) * J * w(igp) * thickness;
fT = (N' * S) * J * w(igp) * thickness;
```

其中 `M` 可理解为 Robin 边界中的矩阵系数，`S` 是边界源项。

### 代码评价

逻辑正确。

建议补充长度检查：

```matlab
if L <= 0
    error('Zero-length boundary edge.');
end
```

但这不是必须改动。

### 考试重点

边界单元最可能考：

```matlab
J = L/2;
N = [(1-s)/2, (1+s)/2];
```

以及 traction 到节点力的积分。

---

## 3.11 `assembleNaturalBC_ThermoMech.m`

### 作用

遍历 `model.natBC` 中所有边界段，把 `Edge2_ThermoMech` 的边界矩阵和边界向量装配到全局矩阵中。

### 核心代码

```matlab
edgeDoFs = [(nodeA-1)*nNodeDoF+1;
            (nodeA-1)*nNodeDoF+2;
            (nodeA-1)*nNodeDoF+3;
            (nodeB-1)*nNodeDoF+1;
            (nodeB-1)*nNodeDoF+2;
            (nodeB-1)*nNodeDoF+3];

bcMatrix(edgeDoFs, edgeDoFs) = bcMatrix(edgeDoFs, edgeDoFs) + edgeMatrix;
bcVector(edgeDoFs) = bcVector(edgeDoFs) + edgeVector;
```

### FEM 意义

这和域内单元装配完全同一个逻辑：

```text
局部边界自由度 -> 全局自由度
局部边界矩阵 -> 全局边界矩阵
局部边界向量 -> 全局边界向量
```

### 代码评价

正确。

注意：`findEdgeOwner` 只检查这条边是否属于某个单元，不严格判断它是不是外边界。因此如果用户把内部边误写进 `natBC`，代码也可能接受。这是建模输入风险，不是当前函数的语法错误。

---

## 3.12 `findEdgeOwner.m`

### 作用

判断一条边 `(nodeA,nodeB)` 是否存在于某个 Tri3 或 Quad4 单元中。

### 核心逻辑

Tri3 边：

```matlab
[1 2], [2 3], [3 1]
```

Quad4 边：

```matlab
[1 2], [2 3], [3 4], [4 1]
```

### 代码评价

用于查找边是正确的。

但它没有区分：

```text
外边界边
内部共享边
```

若要更严格，可以统计同一条边属于几个单元：

```text
只属于 1 个单元 -> 外边界
属于 2 个单元 -> 内部边
```

---

## 3.13 `assembleSystem.m`

### 作用

装配域内所有单元的全局矩阵和全局向量。

### 主流程

```text
1. 遍历 element set
2. 读取 eleName
3. 调用 getElementData
4. 遍历单元
5. 读取节点、材料、载荷
6. 调用单元函数
7. getElementDoFs
8. 装配 globalMatrix 和 globalVector
```

### 核心代码

```matlab
[eleMatrix, eleVector] = eleData.func(eleNodesInfo, material, problem, loads, state);
eleDoFs = getElementDoFs(eleNodeIDs, nNodeDoF);

globalMatrix(eleDoFs, eleDoFs) = globalMatrix(eleDoFs, eleDoFs) + eleMatrix;
globalVector(eleDoFs) = globalVector(eleDoFs) + eleVector;
```

### 代码评价

核心逻辑正确。

建议把空单元集检查提前到 `getElementData` 前面。

原始顺序：

```matlab
eleData = getElementData(eleName, model);
...
if isempty(elesSet)
    continue;
end
```

建议顺序：

```matlab
elesSet = model.elesInfo.(elesSetName{setI});
if isempty(elesSet)
    continue;
end

eleData = getElementData(eleName, model);
```

这样更稳。

### 考试重点

这是最核心的代码之一。考试很可能要求补全装配三行。

---

## 3.14 `buildSystem.m`

### 作用

总装配入口。

先装配域内项：

```matlab
[globalMatrix, globalVector] = assembleSystem(model);
```

再根据 `model.fieldNames` 自动选择自然边界条件函数。

### 热力耦合部分

```matlab
elseif isequal(fieldNames, {'ux','uy','T'})
    [bcMatrix, bcVector] = assembleNaturalBC_ThermoMech(model);
    globalMatrix = globalMatrix + bcMatrix;
    globalVector = globalVector + bcVector;
```

### 代码评价

正确。

---

## 3.15 `buildBCData.m`

### 作用

把 `model.essBC` 转成：

```text
bcData.prescribedDoFs
bcData.prescribedValues
```

### 当前格式

热力耦合中每个节点必须一次性给完整 3 个自由度：

```matlab
model.essBC = {
    [node list], [ux_value, uy_value, T_value]
};
```

### 限制

它不支持只对某个节点施加部分自由度，例如：

```text
只给 T，不给 ux, uy
只给 ux，不给 T
```

所以 `main_mixed_thermomech.m` 才采用了对所有边界节点施加完整精确解的方式。

### 是否需要修改

对于当前自由热膨胀 patch test，不需要修改。

如果以后要做更真实的热力耦合问题，比如：

```text
温度边界很多
机械约束很少
```

则建议扩展为 component-wise BC。

### 考试重点

重复边界条件冲突检查很重要：

```matlab
if isPrescribed(dof)
    if abs(bcData.prescribedValues(dof) - val) > 1e-12
        error('Conflicting essential BC detected at global DoF %d.', dof);
    end
end
```

---

## 3.16 `solveSystem.m`

### 作用

施加本质边界条件并求解：

```text
K u = F
```

### 核心分块

```matlab
activeDoFs = setdiff(allDoFs, prescribedDoFs);
solution(prescribedDoFs) = prescribedValues(prescribedDoFs);

Kaa = globalMatrix(activeDoFs, activeDoFs);
Kap = globalMatrix(activeDoFs, prescribedDoFs);
Fa = globalVector(activeDoFs);
up = solution(prescribedDoFs);

solution(activeDoFs) = Kaa \ (Fa - Kap * up);
reaction = globalMatrix * solution - globalVector;
```

### FEM 意义

将自由度分为：

```text
已知自由度 prescribed
未知自由度 active
```

求解未知部分：

```text
ua = Kaa^{-1}(Fa - Kap*up)
```

### 代码评价

正确。

### 考试重点

这是极高频考点。

---

## 3.17 `postprocessModel.m`

### 作用

遍历所有单元，自动调用对应后处理函数：

```matlab
postFunc = str2func(['post_' eleName]);
```

例如：

```text
Quad4_ThermoMech -> post_Quad4_ThermoMech
Tri3_ThermoMech  -> post_Tri3_ThermoMech
```

### 代码评价

正确。

注意：如果某个单元类型存在单元函数，也必须有对应的 post 函数。

---

## 3.18 `extractNodalFields.m`

### 作用

把全局解向量拆成节点物理量。

对于热力耦合：

```matlab
nodal.ux
nodal.uy
nodal.T
```

### 核心代码

```matlab
for n = 1:numNodes
    base = (n-1)*nNodeDoF;
    for i = 1:nNodeDoF
        nodal.(fieldNames{i})(n) = solution(base+i);
    end
end
```

### 代码评价

正确。

---

## 3.19 `averageFieldToNodes.m`

### 作用

把单元结果平均到节点上，用于绘图。

例如：

```matlab
vmN = averageFieldToNodes(model, results, 'vonMises');
stressN = averageFieldToNodes(model, results, 'stress');
```

### FEM 意义

应力、热流这类量通常是单元量或 Gauss 点量，不是直接的节点未知量。为了画云图，常需要平均到节点。

### 代码评价

正确。

注意：如果混合单元中某个字段维度不一致，会报错。例如 Tri3 的 `u` 是 6 分量，Quad4 的 `u` 是 8 分量，不适合直接平均。

---

## 4. 修正版文件说明

我已生成独立修正版压缩包：

```text
FEM_thermomech_fixed_part3.zip
```

包含：

```text
Tri3_ThermoMech.m
Quad4_ThermoMech.m
post_Quad4_ThermoMech.m
assembleSystem.m
main_mixed_thermomech_fixed.m
README_fixed_part3.txt
```

### 4.1 新增 `Tri3_ThermoMech.m`

用于补齐 mixed mesh 中 Tri3 热力耦合单元。

### 4.2 修正 `Quad4_ThermoMech.m`

修正内容：

```text
1. 使用标准 Q4 节点顺序
2. 加入 Tref 对机械右端项的贡献
3. 加入 detJ <= 0 检查
4. 加入 bodyForce 维度检查
5. 用 J \ (...) 替代 inv(J)*(...)
```

### 4.3 修正 `post_Quad4_ThermoMech.m`

修正内容：

```text
1. 与 Quad4_ThermoMech.m 使用相同 Q4 节点顺序
2. 加入 detJ 检查
3. 用 J \ (...) 替代 inv(J)*(...)
```

### 4.4 修正 `assembleSystem.m`

修正内容：

```text
1. 先跳过空 element set，再解析单元函数
2. 增加 node 和 section 的存在性检查
```

### 4.5 新增 `main_mixed_thermomech_fixed.m`

修正内容：

```text
1. 使用 nTriCols = 4 验证真正 mixed Tri3 + Quad4
2. 设置 Tref = 20 测试 Tref 修正是否一致
3. 绘图函数使用 exist(...) 保护，避免缺少 plot 函数时核心求解无法完成
4. max(...,'all') 改成更兼容的 max(abs(A(:))) 写法
```

---

## 5. 按 FEM 逻辑顺序梳理整套热力耦合代码

## 5.1 建立模型

```matlab
model.fieldNames = {'ux','uy','T'};
```

定义每个节点有三个未知量。

## 5.2 生成网格

主程序生成：

```text
nodesInfo
elesInfo.elesSet1
elesInfo.elesSet2
```

其中：

```text
elesSet1 : Tri3_ThermoMech
elesSet2 : Quad4_ThermoMech
```

## 5.3 设置材料

```matlab
model.section(e,2:6) = [E, nu, alpha, kxx, kyy];
```

其中：

```text
E, nu       -> 机械弹性
alpha       -> 热膨胀系数
kxx, kyy    -> 导热系数
```

## 5.4 设置单元载荷

```matlab
model.eleLoadData{e} = struct('bodyForce', [0;0], 'heatSource', 0.0);
```

包括：

```text
bodyForce  -> 机械体力
heatSource -> 体热源
```

## 5.5 单元计算

每个单元计算：

```text
Kuu, KuT, KTT
fMech, fTherm
```

再嵌入：

```text
9 x 9 Tri3 element matrix
12 x 12 Quad4 element matrix
```

## 5.6 整体装配

```matlab
globalMatrix(eleDoFs, eleDoFs) = globalMatrix(eleDoFs, eleDoFs) + eleMatrix;
globalVector(eleDoFs) = globalVector(eleDoFs) + eleVector;
```

## 5.7 自然边界

如果有 `model.natBC`，则：

```text
Edge2_ThermoMech -> assembleNaturalBC_ThermoMech -> buildSystem
```

## 5.8 本质边界

```text
buildBCData
```

把节点边界条件转成自由度边界条件。

## 5.9 求解

```text
solveSystem
```

分块求解未知自由度。

## 5.10 后处理

```text
extractNodalFields
postprocessModel
averageFieldToNodes
```

得到：

```text
ux, uy, T
strain
thermalStrain
stress
vonMises
gradT
heatFlux
```

---

## 6. 考试高频代码段总结

## 6.1 热力耦合本构

```matlab
strainTh = alpha * (T - Tref) * [1; 1; 0];
stress = D * (strain - strainTh);
```

## 6.2 Q4 Jacobian

```matlab
J = [dNds'*x, dNds'*y;
     dNdt'*x, dNdt'*y];

gradN = J \ [dNds'; dNdt'];
```

## 6.3 Q4 热力耦合矩阵

```matlab
Kuu = Kuu + (Bu' * D * Bu) * detJ * w * thickness;
KTT = KTT + (BT' * Kcond * BT) * detJ * w * thickness;
KuT = KuT - (Bu' * D * epsTh) * N' * detJ * w * thickness;
```

## 6.4 单元自由度分组

```matlab
mechDofs = [1 2 4 5 7 8 10 11];
tempDofs = [3 6 9 12];
```

## 6.5 全局装配

```matlab
eleDoFs = getElementDoFs(eleNodeIDs, nNodeDoF);
globalMatrix(eleDoFs, eleDoFs) = globalMatrix(eleDoFs, eleDoFs) + eleMatrix;
globalVector(eleDoFs) = globalVector(eleDoFs) + eleVector;
```

## 6.6 本质边界求解

```matlab
activeDoFs = setdiff(allDoFs, prescribedDoFs);
solution(activeDoFs) = Kaa \ (Fa - Kap * up);
reaction = globalMatrix * solution - globalVector;
```

## 6.7 自由热膨胀精确解

```matlab
ux = alpha * DeltaT * x;
uy = alpha * DeltaT * y;
T  = Tref + DeltaT;
stress = 0;
```

---

## 7. 复习建议

这批热力耦合代码不要只背文件名，要按下面这条链理解：

```text
T -> thermal strain -> stress -> mechanical equation
```

再对应到代码：

```text
N, BT -> temperature and gradT
Bu -> mechanical strain
D, alpha, Tref -> thermal stress correction
Kuu, KTT, KuT -> coupled element matrix
assembleSystem -> global matrix
solveSystem -> solution
postprocessModel -> stress and heat flux
```

考试时最可能要求你补全：

```text
1. D matrix
2. Bu / BT matrix
3. Kuu / KTT / KuT
4. eleDoFs
5. global assembly
6. prescribed/active DoF solve
7. thermal strain and stress postprocessing
```

---

## 8. 最后判断

如果只运行原始 `main_mixed_thermomech.m` 且保持：

```matlab
nTriCols = 0;
Tref = 0.0;
```

核心求解部分大概率可以通过静态逻辑检查。

但如果用于完整课程代码、混合网格或通用热力耦合问题，则应至少补齐：

```text
Tri3_ThermoMech.m
Quad4_ThermoMech.m 中 Tref 项
Q4 单元与后处理的节点顺序统一
```

因此我建议以后用本次压缩包中的修正版作为热力耦合部分的继续学习版本。
