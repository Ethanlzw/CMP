# FEM MATLAB 脚本第二批代码审查与讲解

## 0. 本批已阅读文件

本批新增 7 个 `.m` 文件：

1. `Tri3_Mech.m`
2. `Quad4_Mech.m`
3. `Tri3_Thermal.m`
4. `Quad4_Thermal.m`
5. `post_Tri3_Thermal.m`
6. `postprocessModel.m`
7. `solveSystem.m`

同时结合第一批中已经阅读的下列文件重新核对接口：

1. `main_mixed_mech.m`
2. `main_mixed_thermal.m`
3. `assembleSystem.m`
4. `buildSystem.m`
5. `buildBCData.m`
6. `getElementData.m`
7. `getElementDoFs.m`
8. `post_Quad4_Mech.m`
9. `post_Quad4_Thermal.m`
10. `getElasticMatrix2D.m`

本报告是第一批报告的续篇，重点补上单元函数、求解器和后处理入口。

## 1. 总体结论

这批文件补齐了最关键的 FEM 计算核心：

```text
Tri3_Mech / Quad4_Mech      负责二维线弹性单元刚度和体力载荷
Tri3_Thermal / Quad4_Thermal 负责二维稳态热传导单元矩阵和热源项
solveSystem                 负责施加本质边界后的线性方程求解
postprocessModel            负责逐单元调用后处理函数
post_Tri3_Thermal           负责三角形热单元的温度梯度和热流后处理
```

纯力学和纯热传导的 FEM 逻辑已经基本闭合。当前代码里有 3 类必须修正的问题：

1. `Quad4_Mech.m` 中 `bodyForce` 未定义，会导致运行错误。
2. `Quad4_Thermal.m` 中 `w` 未定义，会导致运行错误。
3. `Quad4_Thermal.m`、`post_Quad4_Thermal.m`、`post_Quad4_Mech.m` 的 Q4 形函数节点顺序与 `main_mixed_mech.m`、`main_mixed_thermal.m` 中的单元节点顺序不一致，会造成错误的 Jacobian 和错误的梯度结果。

我已生成独立修正版，不覆盖原始脚本：

```text
FEM_fixed_part2.zip
```

压缩包内包含：

```text
Quad4_Mech.m
Quad4_Thermal.m
post_Quad4_Mech.m
post_Quad4_Thermal.m
```

## 2. 必须修正的问题

### 2.1 `Quad4_Mech.m` 中 `bodyForce` 未定义

原始代码第 81 行使用：

```matlab
eleVector = eleVector + Nmat' * bodyForce * detJ * w * thickness;
```

但函数内部没有定义 `bodyForce`。

这会在运行到第一个 Quad4 力学单元时直接报错：

```text
Undefined function or variable 'bodyForce'
```

即使 `main_mixed_mech.m` 中已经写了：

```matlab
model.eleLoadData{e} = struct('bodyForce', [0; 0]);
```

也需要在 `Quad4_Mech.m` 内部读取：

```matlab
bodyForce = [0; 0];
if nargin >= 4 && ~isempty(loads) && isfield(loads,'bodyForce')
    bodyForce = loads.bodyForce(:);
    if numel(bodyForce) ~= 2
        error('Mechanical bodyForce must be a 2-component vector [bx; by].');
    end
end
```

这一处属于必须改动，因为原代码无法完成 Quad4 力学单元载荷积分。

### 2.2 `Quad4_Thermal.m` 中 `w` 未定义

原始代码中：

```matlab
weights = [1; 1; 1; 1];
```

但在 Gauss 循环里没有写：

```matlab
w = weights(igp);
```

后面却使用：

```matlab
eleMatrix = eleMatrix + (B' * Kcond * B) * detJ * w * thickness;
eleVector = eleVector + N * qv * detJ * w * thickness;
```

这也会直接报错。

必须在循环中补上：

```matlab
w = weights(igp);
```

### 2.3 Q4 热单元节点顺序与主程序不一致

主程序生成 Q4 单元时使用：

```matlab
quadElems = [quadElems; eleID, n1, n2, n3, n4];
```

其中：

```text
n1: 左下角
n2: 右下角
n3: 右上角
n4: 左上角
```

这对应标准 Q4 母单元顺序：

```text
1(-1,-1), 2(1,-1), 3(1,1), 4(-1,1)
```

`Quad4_Mech.m` 中的形函数顺序正好符合这个定义：

```matlab
N = 0.25 * [(1 - s) * (1 - t), ...
            (1 + s) * (1 - t), ...
            (1 + s) * (1 + t), ...
            (1 - s) * (1 + t)];
```

但 `Quad4_Thermal.m` 原始代码写成：

```matlab
N = 0.25 * [(1+s)*(1+t);
            (1-s)*(1+t);
            (1-s)*(1-t);
            (1+s)*(1-t)];
```

这个顺序对应：

```text
1(1,1), 2(-1,1), 3(-1,-1), 4(1,-1)
```

它和主程序的实际节点顺序相反。后果包括：

1. `detJ` 可能变成负值。
2. `gradN` 映射错误。
3. 热传导刚度矩阵物理方向错误。
4. 后处理热流方向错误。

必须把 `Quad4_Thermal.m` 改成和 `Quad4_Mech.m` 相同的 Q4 顺序。

### 2.4 `post_Quad4_Mech.m` 与 `post_Quad4_Thermal.m` 也存在 Q4 顺序问题

第一批中已经阅读过这两个文件。当时尚未完整对照单元函数。现在结合本批 `Quad4_Mech.m` 后可以确认：

```matlab
post_Quad4_Mech.m
post_Quad4_Thermal.m
```

使用的形函数导数顺序也和主程序 Q4 顺序不一致。

原始后处理代码在中心点使用：

```matlab
dNds = 0.25 * [(1+t); -(1+t); -(1-t); (1-t)];
dNdt = 0.25 * [(1+s);  (1-s); -(1-s); -(1+s)];
```

该顺序对应的节点排列与主程序不同。因此应力、应变、热流的后处理值会出错。

这属于必须修正，因为即使求解阶段正确，后处理结果也会被错误解释。

### 2.5 建议增加 `detJ` 或面积符号检查

Tri3 和 Quad4 单元都依赖单元节点编号方向。

对于 Tri3：

```matlab
Area = detA / 2;
```

对于 Quad4：

```matlab
detJ = det(J);
```

当前主程序生成的三角形和四边形节点顺序是逆时针，面积和 Jacobian 为正。因此在当前两个 main 文件下可以正常运行。

如果未来导入外部网格，节点编号反向时会导致：

```text
Area < 0
或 detJ < 0
```

建议增加检查：

```matlab
if Area <= 0
    error('Tri3 element has non-positive area. Check node order.');
end
```

```matlab
if detJ <= 0
    error('Quad4 element has non-positive detJ. Check node order or distortion.');
end
```

这属于强烈建议项。当前主程序不一定触发该问题。

## 3. FEM 逻辑顺序总览

补齐本批文件后，完整流程可以写成：

```text
main_mixed_mech / main_mixed_thermal
    1. 定义节点、单元、材料、载荷、边界
    2. buildSystem(model)
        2.1 assembleSystem(model)
            2.1.1 getElementData(eleName, model)
            2.1.2 getElementDoFs(eleNodeIDs, nNodeDoF)
            2.1.3 Tri3_Mech / Quad4_Mech / Tri3_Thermal / Quad4_Thermal
            2.1.4 局部矩阵装配到整体矩阵
        2.2 assembleNaturalBC_Mech / assembleNaturalBC_Thermal
    3. buildBCData(model)
    4. solveSystem(K, F, bcData)
    5. extractNodalFields(model, solution)
    6. postprocessModel(model, solution)
    7. averageFieldToNodes(model, results, fieldName)
    8. plot...
```

核心 FEM 链条：

```text
几何和材料
-> 形函数 N
-> 形函数导数
-> Jacobian
-> B 矩阵
-> 单元矩阵 Ke
-> 整体矩阵 K
-> 本质边界条件
-> 解向量 u 或 T
-> 后处理应变、应力、温度梯度、热流
```

## 4. `Tri3_Mech.m` 审查与讲解

### 4.1 文件作用

`Tri3_Mech.m` 是三节点常应变三角形单元，适用于二维线弹性问题。每个节点有两个自由度：

```text
ux, uy
```

单元自由度顺序为：

```text
[ux1, uy1, ux2, uy2, ux3, uy3]^T
```

输出：

```matlab
eleMatrix  % 6 x 6 单元刚度矩阵
eleVector  % 6 x 1 单元体力等效节点力
info       % 面积、B矩阵、D矩阵等调试信息
```

### 4.2 几何矩阵和面积

代码：

```matlab
A = [1 x1 y1;
     1 x2 y2;
     1 x3 y3];
detA = det(A);
Area = detA / 2;
```

这来自线性三角形形函数的面积坐标推导。对三角形单元：

```text
detA = 2A
```

其中 `A` 是单元面积。

考试重点：三角形节点顺序若反向，`detA` 会变负。当前主程序生成的 Tri3 顺序是逆时针，面积为正。

### 4.3 `b_i` 和 `c_i`

代码：

```matlab
b1 = y2 - y3;  c1 = x3 - x2;
b2 = y3 - y1;  c2 = x1 - x3;
b3 = y1 - y2;  c3 = x2 - x1;
```

它们来自三角形形函数导数：

```text
dNi/dx = bi / detA
dNi/dy = ci / detA
```

这是一类高频手写题。考题可能给三角形三个节点坐标，让你补 `b1,b2,b3,c1,c2,c3`。

### 4.4 B 矩阵

代码：

```matlab
B = [b1 0  b2 0  b3 0;
     0  c1 0  c2 0  c3;
     c1 b1 c2 b2 c3 b3] / detA;
```

对应二维小变形关系：

```text
exx     = du/dx
eyy     = dv/dy
gammaxy = du/dy + dv/dx
```

所以：

```text
strain = B * ue
```

其中：

```text
ue = [ux1 uy1 ux2 uy2 ux3 uy3]^T
```

Tri3 是线性位移单元，因此 `B` 在单元内为常数，应变也为常数。

### 4.5 材料矩阵 D

代码：

```matlab
D = getElasticMatrix2D(E, nu, mechType);
```

`D` 根据 `mechType` 选择平面应力或平面应变。力学单元核心公式为：

```text
stress = D * strain
```

### 4.6 单元刚度矩阵

代码：

```matlab
eleMatrix = thickness * Area * (B' * D * B);
```

对应公式：

```text
Ke = ∫ B^T D B dV
```

对常应变三角形单元，`B` 和 `D` 都为常数，体积积分简化为：

```text
Ke = t * A * B^T D B
```

考试重点：这行非常可能作为填空题出现。

### 4.7 体力一致载荷

代码：

```matlab
eleVector = thickness * Area / 3 * [bodyForce(1);
                                    bodyForce(2);
                                    bodyForce(1);
                                    bodyForce(2);
                                    bodyForce(1);
                                    bodyForce(2)];
```

对应常体力的一致节点力：

```text
fe = ∫ N^T b dV
```

对于 Tri3，三个形函数在单元上的积分相同：

```text
∫ Ni dA = A/3
```

所以每个节点分到相同的体力贡献。

### 4.8 风险点

`Tri3_Mech.m` 中直接写：

```matlab
bodyForce = loads.bodyForce(:);
```

如果某个模型没有提供 `loads.bodyForce`，会报错。当前 `main_mixed_mech.m` 已经为每个元素提供：

```matlab
model.eleLoadData{e} = struct('bodyForce', [0; 0]);
```

所以当前主程序可以通过这一点。为了提高通用性，建议改成默认零体力。

## 5. `Quad4_Mech.m` 审查与讲解

### 5.1 文件作用

`Quad4_Mech.m` 是四节点双线性四边形力学单元。每个节点两个自由度：

```text
ux, uy
```

单元自由度顺序：

```text
[ux1, uy1, ux2, uy2, ux3, uy3, ux4, uy4]^T
```

### 5.2 Q4 形函数

标准节点顺序：

```text
1(-1,-1)
2( 1,-1)
3( 1, 1)
4(-1, 1)
```

代码：

```matlab
N = 0.25 * [(1 - s) * (1 - t), ...
            (1 + s) * (1 - t), ...
            (1 + s) * (1 + t), ...
            (1 - s) * (1 + t)];
```

对应公式：

```text
N1 = 1/4(1-s)(1-t)
N2 = 1/4(1+s)(1-t)
N3 = 1/4(1+s)(1+t)
N4 = 1/4(1-s)(1+t)
```

这与主程序生成的 `[n1,n2,n3,n4]` 一致。

### 5.3 Jacobian

代码：

```matlab
J = [dNds * x, dNds * y;
     dNdt * x, dNdt * y];
detJ = det(J);
```

Jacobian 用于把母单元中的导数映射到实际坐标：

```text
[dN/dx; dN/dy] = J^{-1} [dN/ds; dN/dt]
```

代码中使用：

```matlab
gradN = J \ [dNds; dNdt];
```

这是正确写法，比 `inv(J)*...` 更稳定。

### 5.4 B 矩阵

代码：

```matlab
B = [dNdx(1)      0    dNdx(2)   0      dNdx(3)     0      dNdx(4)    0;
       0      dNdy(1)    0     dNdy(2)    0      dNdy(3)    0       dNdy(4);
     dNdy(1)  dNdx(1)  dNdy(2) dNdx(2)  dNdy(3)  dNdx(3)   dNdy(4)  dNdx(4)];
```

它仍然对应：

```text
exx     = du/dx
eyy     = dv/dy
gammaxy = du/dy + dv/dx
```

Q4 的 `B` 会随 `(s,t)` 改变，因此需要 Gauss 积分。

### 5.5 2 x 2 Gauss 积分

代码使用 4 个积分点：

```matlab
gaussPts = [-1/sqrt(3), -1/sqrt(3);
            -1/sqrt(3),  1/sqrt(3);
             1/sqrt(3), -1/sqrt(3);
             1/sqrt(3),  1/sqrt(3)];
weights = [1; 1; 1; 1];
```

单元刚度：

```matlab
eleMatrix = eleMatrix + (B' * D * B) * detJ * w * thickness;
```

对应公式：

```text
Ke = ∫∫ B^T D B * t * detJ ds dt
```

考试重点：Q4 单元常考 `N`、`dN/ds`、`dN/dt`、`J`、`B`、`2x2 Gauss`。

### 5.6 必须修正点

原始文件缺少：

```matlab
bodyForce = ...
```

所以 `eleVector` 无法计算。修正版已在压缩包中给出。

## 6. `Tri3_Thermal.m` 审查与讲解

### 6.1 文件作用

`Tri3_Thermal.m` 是三节点线性热传导单元。每个节点只有一个自由度：

```text
T
```

单元自由度顺序：

```text
[T1, T2, T3]^T
```

### 6.2 温度梯度矩阵

代码：

```matlab
B = [b1 b2 b3;
     c1 c2 c3] / detA;
```

这里的 `B` 表示温度梯度矩阵：

```text
gradT = [dT/dx; dT/dy] = B * Te
```

其中：

```text
Te = [T1 T2 T3]^T
```

### 6.3 热传导材料矩阵

代码：

```matlab
Kcond = [kxx 0;
         0   kyy];
```

各向同性材料时：

```text
kxx = kyy = k
```

各向异性材料时可以用不同的 `kxx` 和 `kyy`。

### 6.4 单元导热矩阵

代码：

```matlab
eleMatrix = thickness * Area * (B' * Kcond * B);
```

对应弱形式：

```text
Ke = ∫ (gradN)^T k (gradN) dV
```

对 Tri3，温度插值为线性，温度梯度为常数，所以积分简化为面积乘积。

### 6.5 体热源项

代码：

```matlab
eleVector = eleVector + thickness * Area * [1/3; 1/3; 1/3] * qv;
```

对应：

```text
fe = ∫ N^T qv dV
```

对线性三角形：

```text
∫ Ni dA = A/3
```

所以三个节点平均分配。

### 6.6 代码判断

`Tri3_Thermal.m` 的核心公式正确。建议增加面积检查，但在当前主程序网格下没有必要改动。

## 7. `Quad4_Thermal.m` 审查与讲解

### 7.1 文件作用

`Quad4_Thermal.m` 是四节点双线性热传导单元。每个节点一个自由度：

```text
T
```

输出：

```matlab
eleMatrix  % 4 x 4 导热矩阵
eleVector  % 4 x 1 热源等效节点向量
```

### 7.2 正确的 Q4 热单元链条

应采用和 `Quad4_Mech.m` 相同的标准节点顺序：

```text
1(-1,-1), 2(1,-1), 3(1,1), 4(-1,1)
```

正确链条为：

```text
N(s,t)
-> dN/ds, dN/dt
-> J
-> dN/dx, dN/dy
-> B = gradN
-> Ke = ∫ B^T Kcond B * detJ * t ds dt
-> fe = ∫ N^T qv * detJ * t ds dt
```

### 7.3 原始代码问题

问题 1：`w` 未定义。

问题 2：形函数顺序与主程序不一致。

问题 3：使用 `inv(J)`：

```matlab
invJ = inv(J);
gradN = invJ * [dNds'; dNdt'];
```

建议改成：

```matlab
gradN = J \ [dNds; dNdt];
```

这样数值上更稳，也更符合 MATLAB 写法。

### 7.4 修正后核心代码

```matlab
N = 0.25 * [(1-s)*(1-t);
            (1+s)*(1-t);
            (1+s)*(1+t);
            (1-s)*(1+t)];

dNds = 0.25 * [-(1-t),  (1-t),  (1+t), -(1+t)];
dNdt = 0.25 * [-(1-s), -(1+s),  (1+s),  (1-s)];

J = [dNds*x, dNds*y;
     dNdt*x, dNdt*y];

gradN = J \ [dNds; dNdt];
B = gradN;

eleMatrix = eleMatrix + (B' * Kcond * B) * detJ * w * thickness;
eleVector = eleVector + N * qv * detJ * w * thickness;
```

### 7.5 考试重点

热传导 Q4 和力学 Q4 的结构非常相似，区别在于：

```text
力学：B 是 3 x 8，应变 = B * ue
热学：B 是 2 x 4，温度梯度 = B * Te
```

力学单元：

```text
Ke = ∫ B^T D B dV
```

热传导单元：

```text
Ke = ∫ B^T Kcond B dV
```

这两个公式是很高频的概念对比题。

## 8. `solveSystem.m` 审查与讲解

### 8.1 文件作用

`solveSystem.m` 处理线性方程：

```text
K u = F
```

并施加本质边界条件，例如：

```text
ux = 0
uy = 0
T = T0
```

函数接口：

```matlab
[reaction, solution] = solveSystem(globalMatrix, globalVector, bcData)
```

这与两个主程序的调用一致。

### 8.2 自由度划分

代码：

```matlab
totalDoFs = length(globalVector);
allDoFs = 1:totalDoFs;

prescribedDoFs = bcData.prescribedDoFs;
prescribedValues = bcData.prescribedValues;
activeDoFs = setdiff(allDoFs, prescribedDoFs);
```

含义：

```text
prescribedDoFs: 已知自由度
activeDoFs: 未知自由度
```

### 8.3 写入已知值

代码：

```matlab
solution = zeros(totalDoFs,1);
solution(prescribedDoFs) = prescribedValues(prescribedDoFs);
```

注意：`buildBCData.m` 中的 `prescribedValues` 是完整长度向量，所以这里用 `prescribedValues(prescribedDoFs)` 是正确的。

### 8.4 分块求解

代码：

```matlab
Kaa = globalMatrix(activeDoFs, activeDoFs);
Kap = globalMatrix(activeDoFs, prescribedDoFs);
Fa = globalVector(activeDoFs);
up = solution(prescribedDoFs);

solution(activeDoFs) = Kaa \ (Fa - Kap * up);
```

对应分块方程：

```text
[Kaa Kap] [ua] = [Fa]
[Kpa Kpp] [up]   [Fp]
```

其中 `up` 已知，所以：

```text
Kaa ua + Kap up = Fa
ua = Kaa^{-1}(Fa - Kap up)
```

考试重点：这段是边界条件处理最可能考的代码。

### 8.5 反力计算

代码：

```matlab
reaction = globalMatrix * solution - globalVector;
```

含义：完整方程回代后，约束自由度上的残量就是反力或反作用量。

对力学问题：

```text
reaction = 支座反力
```

对热传导问题：

```text
reaction = 定温边界上维持温度所需的等效热流反力
```

### 8.6 代码判断

该函数接口和主程序一致，核心逻辑正确。需要注意输出顺序是：

```text
reaction, solution
```

这个顺序和很多教材或旧代码相反，后续调用时不要写反。

## 9. `postprocessModel.m` 审查与讲解

### 9.1 文件作用

`postprocessModel.m` 是统一后处理入口。它不会自己计算应力或热流，而是根据单元名称自动调用对应后处理函数。

核心代码：

```matlab
postFunc = str2func(['post_' eleName]);
```

例如：

```text
eleName = 'Tri3_Mech'
postFunc = post_Tri3_Mech
```

```text
eleName = 'Quad4_Thermal'
postFunc = post_Quad4_Thermal
```

### 9.2 工作流程

每个单元执行：

```text
1. 读取 eleID
2. 读取 eleNodeIDs
3. 根据节点 ID 找坐标
4. 根据 eleID 找材料
5. getElementDoFs 得到全局自由度
6. 从全局 solution 中取 eleSol
7. 调用 post_xxx 函数
8. 存入 results(resCount)
```

### 9.3 关键代码

```matlab
eleDoFs = getElementDoFs(eleNodeIDs, nNodeDoF);
eleSol = solution(eleDoFs);

results(resCount).value = postFunc(eleNodesInfo, material, problem, eleSol);
```

这段非常重要。它把整体解向量拆回单元局部解，再交给后处理函数。

### 9.4 风险点

风险 1：所有单元集必须有相同每节点自由度数。

代码取：

```matlab
nNodeDoF = elesType{1,3};
```

这和当前框架一致，因为同一个模型内只混合同一物理场的不同单元类型，例如：

```text
Tri3_Mech + Quad4_Mech
Tri3_Thermal + Quad4_Thermal
```

风险 2：`post_` 函数必须存在。

如果模型中有 `Quad4_ThermoMech`，但没有 `post_Quad4_ThermoMech.m`，这里会报错。

风险 3：`model.section` 必须按 `eleID` 提供材料行。

代码：

```matlab
matMask = model.section(:,1) == eleID;
material = model.section(matMask,2:end);
```

所以 `section` 第一列必须是元素 ID。

## 10. `post_Tri3_Thermal.m` 审查与讲解

### 10.1 文件作用

`post_Tri3_Thermal.m` 用于三角形热单元后处理，输出：

```text
centroid   单元中心坐标
T          单元平均温度
gradT      温度梯度
heatFlux   热流向量
```

### 10.2 温度梯度

代码：

```matlab
B = [b1 b2 b3;
     c1 c2 c3] / detA;

gradT = B * eleSol;
```

其中：

```text
eleSol = [T1 T2 T3]^T
```

由于 Tri3 温度场线性，温度梯度在单元内为常数。

### 10.3 Fourier 导热定律

代码：

```matlab
q = -Kcond * gradT;
```

对应公式：

```text
q = -k gradT
```

负号表示热流从高温流向低温。

### 10.4 单元代表温度

代码：

```matlab
out.T = mean(eleSol);
```

对线性三角形，中心点温度等于三个节点温度平均值，因此这个写法合理。

### 10.5 代码判断

核心公式正确。建议增加面积符号检查，但在当前主程序生成的网格下没有必要改动。

## 11. `post_Quad4_Mech.m` 和 `post_Quad4_Thermal.m` 的补充审查

这两个文件属于第一批文件。本批结合 `Quad4_Mech.m` 后重新判断，发现需要修正。

### 11.1 问题来源

主程序和 `Quad4_Mech.m` 使用标准顺序：

```text
1 左下
2 右下
3 右上
4 左上
```

第一批原始后处理文件使用另一套 Q4 形函数顺序。这会导致中心点梯度、应变、应力和热流计算错误。

### 11.2 正确处理方式

后处理必须和单元矩阵使用同一套形函数顺序。

因此修正版 `post_Quad4_Mech.m` 和 `post_Quad4_Thermal.m` 已放入：

```text
FEM_fixed_part2.zip
```

### 11.3 考试重点

求解和后处理必须使用同一套节点顺序。常见错误是：

```text
Ke 用一种节点顺序
post 用另一种节点顺序
```

这样位移结果看似能求出，但应力或热流会错。

## 12. 关键代码段与考试重点汇总

### 12.1 Tri3 力学 B 矩阵

```matlab
B = [b1 0  b2 0  b3 0;
     0  c1 0  c2 0  c3;
     c1 b1 c2 b2 c3 b3] / detA;
```

必须理解：

```text
strain = B * ue
```

### 12.2 Tri3 热传导梯度矩阵

```matlab
B = [b1 b2 b3;
     c1 c2 c3] / detA;
```

必须理解：

```text
gradT = B * Te
```

### 12.3 力学刚度矩阵

```matlab
Ke = thickness * Area * (B' * D * B);
```

或 Q4 积分形式：

```matlab
eleMatrix = eleMatrix + (B' * D * B) * detJ * w * thickness;
```

### 12.4 热传导矩阵

```matlab
Ke = thickness * Area * (B' * Kcond * B);
```

或 Q4 积分形式：

```matlab
eleMatrix = eleMatrix + (B' * Kcond * B) * detJ * w * thickness;
```

### 12.5 Q4 Jacobian 和导数映射

```matlab
J = [dNds*x, dNds*y;
     dNdt*x, dNdt*y];

gradN = J \ [dNds; dNdt];
```

必须理解：

```text
自然坐标导数 -> 实际坐标导数
```

### 12.6 整体自由度映射

```matlab
eleDoFs = getElementDoFs(eleNodeIDs, nNodeDoF);
```

对应公式：

```matlab
globalDof = (nodeID - 1) * nNodeDoF + localDof;
```

### 12.7 本质边界条件分块求解

```matlab
solution(activeDoFs) = Kaa \ (Fa - Kap * up);
```

必须掌握分块矩阵推导。

### 12.8 反力

```matlab
reaction = globalMatrix * solution - globalVector;
```

这是支座反力或边界反作用量的来源。

## 13. 按文件列出最终判断

| 文件 | 判断 | 是否必须改 |
|---|---:|---:|
| `Tri3_Mech.m` | 公式正确，建议增加默认体力和面积检查 | 否 |
| `Quad4_Mech.m` | `bodyForce` 未定义 | 是 |
| `Tri3_Thermal.m` | 公式正确，建议增加面积检查 | 否 |
| `Quad4_Thermal.m` | `w` 未定义，Q4 顺序错误 | 是 |
| `solveSystem.m` | 接口和主程序一致，逻辑正确 | 否 |
| `postprocessModel.m` | 通用后处理入口正确 | 否 |
| `post_Tri3_Thermal.m` | 公式正确 | 否 |
| `post_Quad4_Mech.m` | Q4 顺序与主程序不一致 | 是 |
| `post_Quad4_Thermal.m` | Q4 顺序与主程序不一致 | 是 |

## 14. 第二批修正版文件说明

修正版压缩包：

```text
FEM_fixed_part2.zip
```

包含 4 个文件：

```text
Quad4_Mech.m
Quad4_Thermal.m
post_Quad4_Mech.m
post_Quad4_Thermal.m
```

使用方法：

1. 备份原始脚本。
2. 解压 `FEM_fixed_part2.zip`。
3. 将其中 4 个 `.m` 文件复制到 MATLAB 工作目录。
4. 覆盖原文件前建议先保存原版副本。
5. 重新运行 `main_mixed_mech.m` 或 `main_mixed_thermal_fixed.m`。

## 15. 完整知识链条总结

纯力学脚本对应的 FEM 链条：

```text
节点位移 ue
-> B 矩阵
-> strain = B ue
-> stress = D strain
-> Ke = ∫ B^T D B dV
-> 组装 Ku = F
-> 施加位移边界
-> 求解 u
-> 回代 reaction
-> 后处理 stress 和 von Mises
```

纯热传导脚本对应的 FEM 链条：

```text
节点温度 Te
-> gradT = B Te
-> q = -Kcond gradT
-> Ke = ∫ B^T Kcond B dV
-> 组装 KT = Q
-> 施加定温边界
-> 求解 T
-> 回代热反力
-> 后处理温度梯度和热流
```

这两条链条结构高度相似。考试复习时可以统一记成：

```text
场变量插值
-> 对场变量求导
-> 材料关系
-> 弱形式积分
-> 单元矩阵
-> 整体装配
-> 边界条件
-> 求解
-> 后处理
```

## 16. 下一批建议上传的文件

目前仍未看到绘图函数：

```text
plotMesh.m
plotNodalScalar.m
plotDeformedMeshContour.m
plotElementScalar.m
plotDeformedMesh.m
```

如果后续还有热力耦合或压电部分，也建议继续上传：

```text
Tri3_ThermoMech.m
Quad4_ThermoMech.m
post_Tri3_ThermoMech.m
post_Quad4_ThermoMech.m
assembleNaturalBC_ThermoMech.m
Edge2_ThermoMech.m
```

这样可以继续检查热力耦合矩阵分块、热膨胀等效载荷、温度场和力学场的接口一致性。
