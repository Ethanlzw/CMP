# CMP FEM Exam Summary

## 1. 考试总逻辑

FEM 解题主线：
```math
\text{physical problem}
\rightarrow
\text{weak form}
\rightarrow
u \approx Nu_e
\rightarrow
\varepsilon = Bu_e
\rightarrow
\sigma = D\varepsilon
\rightarrow
K^e u^e = F^e
\rightarrow
Ku = F
\rightarrow
\text{BC}
\rightarrow
\text{solution}
\rightarrow
\text{postprocess}
```

对热力耦合问题，未知量是 displacement 和 temperature，边界也分成位移边界、力边界、温度边界、热流边界。

考试里最可能考 4 类：

1. 给一个 Tri3 或 Quad4，推导或补全 element stiffness matrix / force vector。
2. 给代码框架，补全 `Tri3_Mech`、`Quad4_Mech`、`Tri3_Thermal`、`Quad4_ThermoMech` 等关键段。
3. 修改 BC 或 loading 后重新模拟。
4. 在 COMSOL 中做 2D 或 3D thermomechanical case，并和 MATLAB 结果对比。

---

## 2. 必背公式总表

### 2.1 机械场公式

二维小变形应变：

```math
\varepsilon =
\begin{bmatrix}
\varepsilon_{xx}\\
\varepsilon_{yy}\\
\gamma_{xy}
\end{bmatrix}
=
\begin{bmatrix}
u_{,x}\\
v_{,y}\\
u_{,y}+v_{,x}
\end{bmatrix}
```

注意这里用的是 engineering shear strain：

```math
\gamma_{xy}=2\varepsilon_{xy}
```

材料关系：

```math
\sigma = D\varepsilon
```

平面应力：

```math
D =
\frac{E}{1-\nu^2}
\begin{bmatrix}
1 & \nu & 0\\
\nu & 1 & 0\\
0 & 0 & \frac{1-\nu}{2}
\end{bmatrix}
```

平面应变：

```math
D =
\frac{E}{(1+\nu)(1-2\nu)}
\begin{bmatrix}
1-\nu & \nu & 0\\
\nu & 1-\nu & 0\\
0 & 0 & \frac{1-2\nu}{2}
\end{bmatrix}
```

机械单元刚度：

```math
K_{uu}^e=\int_{\Omega_e} B_u^T D B_u \, t\, dA
```

机械体力向量：

```math
F_u^e=\int_{\Omega_e} N_u^T b \, t\, dA
```

边界 traction 向量：

```math
F_t^e=\int_{\Gamma_t^e} N_u^T \bar t \, t\, d\Gamma
```

### 2.2 热传导公式

Fourier law：

```math
q=-K\nabla T
```

稳态热传导单元矩阵：

```math
K_{TT}^e=\int_{\Omega_e} B_T^T K B_T \, t\, dA
```

体热源向量：

```math
F_T^e=\int_{\Omega_e} N_T^T r \, t\, dA
```

热流边界向量：

```math
F_q^e=\int_{\Gamma_q^e} N_T^T \bar q_n \, t\, d\Gamma
```

Robin 边界：

```math
K_\Gamma^e=\int_{\Gamma_q^e} N_T^T hN_T \,t\,d\Gamma
```

```math
F_\Gamma^e=\int_{\Gamma_q^e} N_T^T hT_\infty \,t\,d\Gamma
```

### 2.3 热力耦合公式

热应变：

```math
\varepsilon_{th}
=
\alpha (T-T_{ref})
\begin{bmatrix}
1\\
1\\
0
\end{bmatrix}
```

热弹性本构：

```math
\sigma = D(\varepsilon-\varepsilon_{th})
```

单元热力耦合矩阵：

```math
K_{uT}^e
=
-\int_{\Omega_e}
B_u^T D
\left(
\alpha
\begin{bmatrix}
1\\1\\0
\end{bmatrix}
N_T
\right)
t\,dA
```

热力耦合单元系统：

```math
\begin{bmatrix}
K_{uu}^e & K_{uT}^e\\
0 & K_{TT}^e
\end{bmatrix}
\begin{bmatrix}
d_e\\
\theta_e
\end{bmatrix}
=
\begin{bmatrix}
f_u^e\\
f_T^e
\end{bmatrix}
```

整体系统：

```math
\begin{bmatrix}
K_{uu} & K_{uT}\\
0 & K_{TT}
\end{bmatrix}
\begin{bmatrix}
d\\
\theta
\end{bmatrix}
=
\begin{bmatrix}
f_u\\
f_T
\end{bmatrix}
```

含义：先由 $`K_{TT}\theta=f_T`$ 求温度，再由温度场产生等效热膨胀载荷，进入机械方程。

---

## 3. Tri3 单元考试公式

### 3.1 Tri3 形函数

三角形节点：

```math
(x_1,y_1),\quad (x_2,y_2),\quad (x_3,y_3)
```

面积：

```math
A=
\frac12
\det
\begin{bmatrix}
1&x_1&y_1\\
1&x_2&y_2\\
1&x_3&y_3
\end{bmatrix}
```

常用定义：

```math
b_1=y_2-y_3,\quad b_2=y_3-y_1,\quad b_3=y_1-y_2
```

```math
c_1=x_3-x_2,\quad c_2=x_1-x_3,\quad c_3=x_2-x_1
```

形函数导数：

```math
N_{i,x}=\frac{b_i}{2A}
```

```math
N_{i,y}=\frac{c_i}{2A}
```

Tri3 是常应变单元，所以 $`B`$ 在单元内部是常数。

### 3.2 Tri3 机械 $`B_u`$

```math
B_u=
\frac{1}{2A}
\begin{bmatrix}
b_1&0&b_2&0&b_3&0\\
0&c_1&0&c_2&0&c_3\\
c_1&b_1&c_2&b_2&c_3&b_3
\end{bmatrix}
```

机械刚度：

```math
K_{uu}^e=B_u^TDB_u\,tA
```

体力：

```math
F_u^e=
\frac{tA}{3}
\begin{bmatrix}
b_x\\b_y\\b_x\\b_y\\b_x\\b_y
\end{bmatrix}
```

### 3.3 Tri3 热学 $`B_T`$

```math
B_T=
\frac{1}{2A}
\begin{bmatrix}
b_1&b_2&b_3\\
c_1&c_2&c_3
\end{bmatrix}
```

热刚度：

```math
K_{TT}^e=B_T^T K B_T\,tA
```

体热源：

```math
F_T^e=
\frac{r\,tA}{3}
\begin{bmatrix}
1\\1\\1
\end{bmatrix}
```

### 3.4 Tri3 热力耦合

如果温度也用同一个 Tri3 形函数：

```math
T=N_T\theta_e
```

其中：

```math
N_T=\begin{bmatrix}N_1&N_2&N_3\end{bmatrix}
```

则：

```math
K_{uT}^e
=
-\int_{\Omega_e}
B_u^TD
\left(
\alpha
\begin{bmatrix}
1\\1\\0
\end{bmatrix}
N_T
\right)t\,dA
```

Tri3 里 $`B_u`$ 常数，若用 centroid 积分：

```math
N_1=N_2=N_3=\frac13
```

```math
K_{uT}^e
\approx
-B_u^T D
\left(
\alpha
\begin{bmatrix}
1\\1\\0
\end{bmatrix}
\begin{bmatrix}
\frac13&\frac13&\frac13
\end{bmatrix}
\right)tA
```

---

## 4. Q4 单元考试公式

### 4.1 Q4 形函数

自然坐标：

```math
(s,t)\in[-1,1]\times[-1,1]
```

标准节点顺序：

1. $`(-1,-1)`$
2. $`(1,-1)`$
3. $`(1,1)`$
4. $`(-1,1)`$

形函数：

```math
N_1=\frac14(1-s)(1-t)
```

```math
N_2=\frac14(1+s)(1-t)
```

```math
N_3=\frac14(1+s)(1+t)
```

```math
N_4=\frac14(1-s)(1+t)
```

### 4.2 Q4 Jacobian

几何映射：

```math
x=\sum_{i=1}^4N_ix_i
```

```math
y=\sum_{i=1}^4N_iy_i
```

代码中若定义：

```matlab
J = [dNds'; dNdt'] * coords;
```

其中：

```matlab
coords = [x1 y1;
          x2 y2;
          x3 y3;
          x4 y4];
```

则：

```math
J=
\begin{bmatrix}
x_{,s} & y_{,s}\\
x_{,t} & y_{,t}
\end{bmatrix}
```

导数变换：

```math
\begin{bmatrix}
N_{,x}\\
N_{,y}
\end{bmatrix}
=
J^{-1}
\begin{bmatrix}
N_{,s}\\
N_{,t}
\end{bmatrix}
```

代码中应写：

```matlab
gradN = J \ [dNds'; dNdt'];
dNdx = gradN(1,:)';
dNdy = gradN(2,:)';
```

`J \ rhs` 和 `inv(J)*rhs` 数学上等价，但代码里推荐 `J \ rhs`，数值更稳定。

### 4.3 Q4 机械 $`B_u`$

```math
B_u=
\begin{bmatrix}
N_{1,x}&0&N_{2,x}&0&N_{3,x}&0&N_{4,x}&0\\
0&N_{1,y}&0&N_{2,y}&0&N_{3,y}&0&N_{4,y}\\
N_{1,y}&N_{1,x}&N_{2,y}&N_{2,x}&N_{3,y}&N_{3,x}&N_{4,y}&N_{4,x}
\end{bmatrix}
```

Q4 热梯度矩阵：

```math
B_T=
\begin{bmatrix}
N_{1,x}&N_{2,x}&N_{3,x}&N_{4,x}\\
N_{1,y}&N_{2,y}&N_{3,y}&N_{4,y}
\end{bmatrix}
```

### 4.4 Q4 高斯积分

标准 $`2\times2`$ Gauss 点：

```math
s,t=\pm\frac1{\sqrt3}
```

权重：

```math
w=1
```

Q4 机械刚度：

```math
K_{uu}^e=
\sum_{g=1}^4
B_{u,g}^TDB_{u,g}\det(J_g)w_g t
```

Q4 热刚度：

```math
K_{TT}^e=
\sum_{g=1}^4
B_{T,g}^T K B_{T,g}\det(J_g)w_g t
```

Q4 热力耦合：

```math
K_{uT}^e
=
-\sum_{g=1}^4
B_{u,g}^TD
\left(
\alpha
\begin{bmatrix}
1\\1\\0
\end{bmatrix}
N_{T,g}
\right)
\det(J_g)w_g t
```

因为 Q4 的导数在单元内变化，所以需要数值积分。

---

## 5. 对应脚本清单

### 5.1 主脚本

| 脚本 | 作用 |
|---|---|
| `main_mixed_mech.m` | 纯力学，Tri3 + Quad4 |
| `main_mixed_thermal.m` | 纯热传导，Tri3 + Quad4 |
| `main_mixed_thermomech.m` | 热力耦合 |
| `main_validate_mech_patch.m` | 机械 patch test |
| `main_validate_thermal_fourier.m` | 热传导验证 |
| `main_validate_thermomech_free_expansion.m` | 自由热膨胀验证 |

### 5.2 核心框架脚本

| 脚本 | 作用 |
|---|---|
| `getDofInfo.m` | 根据 `model.fieldNames` 决定自由度顺序 |
| `getElementData.m` | 注册元素类型，返回元素函数句柄 |
| `getElementDof.m` / `getElementDoFs.m` | 节点号转整体自由度号 |
| `assembleSystem.m` | 循环单元，组装整体矩阵 |
| `buildSystem.m` | 组装域单元和自然边界 |
| `buildBCData.m` | 本质边界转 known / unknown DOF |
| `solveSystem.m` | 分块求解 |
| `postprocessModel.m` | 后处理 |
| `extractNodalFields.m` | 从总解向量提取 `ux,uy,T` |
| `plotMesh.m` / `plotNodalScalar.m` / `plotDeformedMesh.m` | 绘图 |

### 5.3 单元脚本

| 脚本 | 主要输出 |
|---|---|
| `Tri3_Mech.m` | `Ke`, `Fe`，尺寸 $`6\times6`$, $`6\times1`$ |
| `Quad4_Mech.m` | `Ke`, `Fe`，尺寸 $`8\times8`$, $`8\times1`$ |
| `Tri3_Thermal.m` | `Ke`, `Fe`，尺寸 $`3\times3`$, $`3\times1`$ |
| `Quad4_Thermal.m` | `Ke`, `Fe`，尺寸 $`4\times4`$, $`4\times1`$ |
| `Tri3_ThermoMech.m` | block matrix，尺寸 $`9\times9`$ |
| `Quad4_ThermoMech.m` | block matrix，尺寸 $`12\times12`$ |

---

## 6. 代码块范例

### 6.1 `assembleSystem.m` 核心循环

```matlab
globalMatrix = zeros(nodesNum*nNodeDoF);
globalVector = zeros(nodesNum*nNodeDoF,1);

for setI = 1:elesSetNum

    elesSet = model.elesInfo.(model.elesSetName{setI});

    if isempty(elesSet)
        continue;
    end

    eleName = model.elesType{setI,2};
    eleData = getElementData(eleName, model);

    for e = 1:size(elesSet,1)

        eleID = elesSet(e,1);
        eleNodes = elesSet(e,2:end);

        coords = model.nodesInfo(eleNodes,2:3);

        [Ke, Fe] = eleData.func(coords, model, eleID);

        edofs = getElementDof(eleNodes, nNodeDoF);

        globalMatrix(edofs,edofs) = globalMatrix(edofs,edofs) + Ke;
        globalVector(edofs) = globalVector(edofs) + Fe;
    end
end
```

考试补全点通常是：

```matlab
globalMatrix(edofs,edofs) = globalMatrix(edofs,edofs) + Ke;
globalVector(edofs) = globalVector(edofs) + Fe;
```

### 6.2 `solveSystem.m` 标准分块

```matlab
function [u, F] = solveSystem(K, u_bc, Fknown)

n = size(K,1);
u_bc = u_bc(:);
Fknown = Fknown(:);

known_u = find(~isnan(u_bc));
un_u    = find(isnan(u_bc));

u = nan(n,1);
u(known_u) = u_bc(known_u);

K_uu = K(un_u, un_u);
K_uk = K(un_u, known_u);

F_u = Fknown(un_u);
u_k = u(known_u);

u(un_u) = K_uu \ (F_u - K_uk*u_k);

F = K*u - Fknown;
end
```

注意反力更推荐写成：

```matlab
reaction = K*u - Fknown;
```

因为 `K*u` 是内力平衡项，`Fknown` 是外载荷项。

### 6.3 `Quad4_Mech.m` 核心代码

```matlab
function [Ke, Fe] = Quad4_Mech(coords, model, eleID)

E  = model.section(eleID,2);
nu = model.section(eleID,3);
t  = model.problem.thickness;

if strcmp(model.problem.mechType,'planeStress')
    D = E/(1-nu^2) * [1 nu 0;
                      nu 1 0;
                      0 0 (1-nu)/2];
else
    D = E/((1+nu)*(1-2*nu)) * [1-nu nu 0;
                               nu 1-nu 0;
                               0 0 (1-2*nu)/2];
end

gauss = [-1/sqrt(3), 1/sqrt(3)];
Ke = zeros(8,8);
Fe = zeros(8,1);

bodyForce = [0;0];
if isfield(model,'eleLoadData')
    bodyForce = model.eleLoadData{eleID}.bodyForce(:);
end

for i = 1:2
    for j = 1:2

        s = gauss(i);
        r = gauss(j);

        N = 1/4 * [(1-s)*(1-r);
                   (1+s)*(1-r);
                   (1+s)*(1+r);
                   (1-s)*(1+r)];

        dNds = 1/4 * [-(1-r);
                       (1-r);
                       (1+r);
                      -(1+r)];

        dNdr = 1/4 * [-(1-s);
                      -(1+s);
                       (1+s);
                       (1-s)];

        J = [dNds'; dNdr'] * coords;
        detJ = det(J);

        if detJ <= 0
            error('Invalid Quad4 element: detJ <= 0. Check node order.');
        end

        gradN = J \ [dNds'; dNdr'];
        dNdx = gradN(1,:)';
        dNdy = gradN(2,:)';

        B = zeros(3,8);
        for a = 1:4
            col = 2*a-1;
            B(:,col:col+1) = [dNdx(a) 0;
                              0 dNdy(a);
                              dNdy(a) dNdx(a)];
        end

        Nu = zeros(2,8);
        for a = 1:4
            Nu(:,2*a-1:2*a) = N(a)*eye(2);
        end

        Ke = Ke + B' * D * B * detJ * t;
        Fe = Fe + Nu' * bodyForce * detJ * t;
    end
end
end
```

### 6.4 `Tri3_Mech.m` 核心代码

```matlab
function [Ke, Fe] = Tri3_Mech(coords, model, eleID)

E  = model.section(eleID,2);
nu = model.section(eleID,3);
t  = model.problem.thickness;

if strcmp(model.problem.mechType,'planeStress')
    D = E/(1-nu^2) * [1 nu 0;
                      nu 1 0;
                      0 0 (1-nu)/2];
else
    D = E/((1+nu)*(1-2*nu)) * [1-nu nu 0;
                               nu 1-nu 0;
                               0 0 (1-2*nu)/2];
end

x = coords(:,1);
y = coords(:,2);

A = 0.5 * det([1 x(1) y(1);
               1 x(2) y(2);
               1 x(3) y(3)]);

if A <= 0
    error('Invalid Tri3 element: area <= 0. Check node order.');
end

b = [y(2)-y(3);
     y(3)-y(1);
     y(1)-y(2)];

c = [x(3)-x(2);
     x(1)-x(3);
     x(2)-x(1)];

B = 1/(2*A) * [b(1) 0 b(2) 0 b(3) 0;
               0 c(1) 0 c(2) 0 c(3);
               c(1) b(1) c(2) b(2) c(3) b(3)];

Ke = B' * D * B * t * A;

bodyForce = [0;0];
if isfield(model,'eleLoadData')
    bodyForce = model.eleLoadData{eleID}.bodyForce(:);
end

Fe = t*A/3 * [bodyForce;
              bodyForce;
              bodyForce];
end
```

### 6.5 `Quad4_ThermoMech.m` 核心结构

```matlab
Ke = zeros(12,12);
Fe = zeros(12,1);

% DOF order per node:
% [ux uy T]
% element vector:
% [ux1 uy1 T1 ux2 uy2 T2 ux3 uy3 T3 ux4 uy4 T4]

for each Gauss point

    % 1. shape functions N
    % 2. dNds, dNdt
    % 3. J, detJ
    % 4. dNdx, dNdy
    % 5. Bu, BT, NT

    Kuu = Bu' * D * Bu * detJ * t;
    KTT = BT' * Kcond * BT * detJ * t;

    m = [1;1;0];
    KuT = -Bu' * D * (alpha * m * NT) * detJ * t;

    % local block dof indices
    uDofs = [1 2 4 5 7 8 10 11];
    tDofs = [3 6 9 12];

    Ke(uDofs,uDofs) = Ke(uDofs,uDofs) + Kuu;
    Ke(tDofs,tDofs) = Ke(tDofs,tDofs) + KTT;
    Ke(uDofs,tDofs) = Ke(uDofs,tDofs) + KuT;

end
```

如果 `fieldNames={'ux','uy','T'}`，每个节点 DOF 顺序必须是：

```matlab
ux, uy, T
```

因此 4 节点 Q4 热力耦合总 DOF 顺序是：

```matlab
[ux1 uy1 T1 ux2 uy2 T2 ux3 uy3 T3 ux4 uy4 T4]
```

---

## 7. 解题步骤模板

### 7.1 题型 A：推导某个单元刚度矩阵

标准答案结构：

1. 写主变量插值：

```math
u=N_u d_e
```

或：

```math
T=N_T\theta_e
```

2. 求导得到：

```math
\varepsilon=B_ud_e
```

或：

```math
\nabla T=B_T\theta_e
```

3. 代入材料关系：

```math
\sigma=D\varepsilon
```

```math
q=-K\nabla T
```

4. 代入弱形式。

5. 得到：

```math
K_{uu}^e=\int B_u^TDB_u\,t\,dA
```

```math
K_{TT}^e=\int B_T^TKB_T\,t\,dA
```

6. 若是 Q4，转换到自然坐标并用 Gauss 积分：

```math
dA=\det(J)\,dsdt
```

7. 写最终矩阵或代码求和形式。

### 7.2 题型 B：补全代码

先判断题目给的是哪一层。

#### 单元层

重点补：

```matlab
N
dNds, dNdt
J
detJ
gradN
B
Ke = Ke + ...
Fe = Fe + ...
```

#### 装配层

重点补：

```matlab
edofs = getElementDof(eleNodes, nNodeDoF);
globalMatrix(edofs,edofs) = globalMatrix(edofs,edofs) + Ke;
globalVector(edofs) = globalVector(edofs) + Fe;
```

#### 求解层

重点补：

```matlab
known = find(~isnan(u_bc));
unknown = find(isnan(u_bc));
u(unknown) = K(unknown,unknown) \ ...
```

#### 后处理层

重点补：

```matlab
strain = B * ue;
stress = D * strain;
gradT  = BT * Te;
heatFlux = -K * gradT;
```

### 7.3 题型 C：测试 element

#### Mechanical patch test

给线性位移场：

```math
u=a_0+a_1x+a_2y
```

```math
v=b_0+b_1x+b_2y
```

应变应为常数。Tri3 和 Q4 都应能给出常应变。

#### Thermal linear patch test

给线性温度场：

```math
T=a+bx+cy
```

温度梯度：

```math
\nabla T=
\begin{bmatrix}
b\\c
\end{bmatrix}
```

热流：

```math
q=-K\nabla T
```

应为常数。

#### Free thermal expansion

若没有机械约束，均匀升温：

```math
T-T_{ref}=\Delta T
```

自由热膨胀时应力应为零：

```math
\sigma=0
```

这是 thermomechanical 最重要验证之一。

#### Restrained thermal expansion

若完全约束位移：

```math
\varepsilon=0
```

则：

```math
\sigma=-D\varepsilon_{th}
```

平面应力下：

```math
\sigma_{xx}=\sigma_{yy}
=
-\frac{E\alpha\Delta T}{1-\nu}
```

---

## 8. 修改 BC 或 loading 的代码模式

### 8.1 纯力学，固定左边界

```matlab
leftNodes = find(abs(nodes(:,2)-0) < tol);

model.essBC = {
    leftNodes, [0,0]
};
```

### 8.2 右边界施加集中力

```matlab
rightNodes = find(abs(nodes(:,2)-L) < tol);

model.nodalVector = zeros(numNodes*2,1);

P = -1000;
for i = 1:numel(rightNodes)
    node = rightNodes(i);
    dofUy = (node-1)*2 + 2;
    model.nodalVector(dofUy) = model.nodalVector(dofUy) + P/numel(rightNodes);
end
```

### 8.3 热传导，左右定温

```matlab
leftNodes  = find(abs(nodes(:,2)-0) < tol);
rightNodes = find(abs(nodes(:,2)-L) < tol);

model.essBC = {
    leftNodes,  [100];
    rightNodes, [0]
};
```

### 8.4 热力耦合，左边固定并给温度

若 `fieldNames={'ux','uy','T'}`：

```matlab
model.essBC = {
    leftNodes, [0,0,100];
    rightNodes, [nan,nan,0]
};
```

更稳妥做法是分别处理每个 field，避免 `nan` 被 `buildBCData` 误处理。考试中写概念时可以说：只约束 `T`，不约束 `ux,uy`。

---

## 9. Pei COMSOL thermomechanical 2D case

### 9.1 COMSOL 建模步骤

1. Model Wizard。
2. Space Dimension 选 **2D**。
3. Physics 选：
   1. Solid Mechanics
   2. Heat Transfer in Solids
   3. Multiphysics 里添加 Thermal Expansion
4. Study 选 Stationary。
5. Geometry 建矩形。
6. Materials：
   1. $`E`$
   2. $`\nu`$
   3. $`\alpha`$
   4. $`k`$
   5. density 如果 stationary 可不作为主参数
7. Solid Mechanics：
   1. 左边固定或指定 displacement
   2. 右边施加 boundary load / traction
8. Heat Transfer：
   1. 左边温度 $`T_{hot}`$
   2. 右边温度 $`T_{cold}`$
   3. 其他边默认 insulation
9. Multiphysics：
   1. Thermal Expansion 中设置 $`T_{ref}`$
10. Mesh：
    1. 用 Free Triangular 或 Mapped Quad
    2. 若要和 MATLAB Q4 对比，尽量用 structured quadrilateral mesh
11. Compute。
12. Results：
    1. Temperature
    2. Displacement magnitude
    3. Stress, von Mises
    4. reaction force

### 9.2 和 MATLAB 对比时必须统一

要统一：

```math
E,\nu,\alpha,k,t,T_{ref}
```

还要统一：

1. plane stress 或 plane strain
2. 几何尺寸
3. 热边界
4. 位移约束
5. 载荷方向
6. 网格密度
7. 结果取点位置

2D COMSOL 默认不等同于 MATLAB plane stress，需要检查 Solid Mechanics 的设置。若 MATLAB 用 plane stress，COMSOL 里要选对应 2D formulation。

---

## 10. Pei COMSOL thermomechanical 3D case

### 10.1 建模步骤

1. Space Dimension 选 **3D**。
2. Geometry 建 Block。
3. Physics：
   1. Solid Mechanics
   2. Heat Transfer in Solids
   3. Thermal Expansion
4. Material 输入：
   1. $`E`$
   2. $`\nu`$
   3. $`\alpha`$
   4. $`k`$
5. Heat Transfer：
   1. 一侧给 $`T_{hot}`$
   2. 另一侧给 $`T_{cold}`$
   3. 其他面 insulation
6. Solid Mechanics：
   1. 一端 fixed constraint
   2. 另一端可 free 或施加载荷
7. Mesh：
   1. Free Tetrahedral 或 Swept Hexahedral
   2. 若想类似 Q4/Q8 逻辑，优先 Swept mesh
8. Study：
   1. Stationary
9. Results：
   1. 3D displacement
   2. stress
   3. temperature
   4. thermal strain

### 10.2 3D 与 2D 的主要差异

3D 中应变有 6 个分量：

```math
\varepsilon =
[\varepsilon_{xx},\varepsilon_{yy},\varepsilon_{zz},
\gamma_{xy},\gamma_{yz},\gamma_{xz}]^T
```

热应变：

```math
\varepsilon_{th}
=
\alpha(T-T_{ref})
[1,1,1,0,0,0]^T
```

所以 3D 不能直接拿 2D 的 $`D`$ 和 $`[1,1,0]^T`$ 套进去。

---

## 11. 易错点清单

### 11.1 节点顺序

Tri3 必须保持逆时针。若面积：

```math
A\le 0
```

说明节点顺序反了。

Q4 标准顺序：

```math
1\rightarrow2\rightarrow3\rightarrow4
```

必须绕单元一圈，通常逆时针。若：

```math
\det J \le 0
```

说明节点顺序或单元形状有问题。

不要用：

```matlab
detJ = abs(det(J));
```

这会掩盖错误的节点顺序。

### 11.2 Q4 导数变换

如果：

```matlab
J = [dNds'; dNdt'] * coords;
```

则：

```matlab
gradN = J \ [dNds'; dNdt'];
```

不要写成：

```matlab
gradN = [dNds'; dNdt'] / J;
```

也不要在没有确认 $`J`$ 定义的情况下随便转置。

### 11.3 `model.fieldNames` 顺序

`model.fieldNames` 是自由度顺序的核心。

常见顺序：

```matlab
{'T'}
{'ux','uy'}
{'ux','uy','T'}
```

如果 `fieldNames={'ux','uy','T'}`，每个节点就是 3 个 DOF。若在单元矩阵里按 `[T ux uy]` 写，装配一定错。

### 11.4 混合单元必须共享同一个 DOF layout

一个模型中可以混合：

```matlab
Tri3_Mech + Quad4_Mech
```

也可以混合：

```matlab
Tri3_Thermal + Quad4_Thermal
```

不能直接混合：

```matlab
Tri3_Mech + Quad4_Thermal
```

除非框架专门支持 field-superset。

### 11.5 本质边界不足

机械问题如果没有足够位移约束，会出现 rigid body mode，矩阵奇异。

2D 弹性至少要消除：

1. $`x`$ 平移
2. $`y`$ 平移
3. 平面转动

典型错误：

```matlab
Warning: Matrix is singular
```

### 11.6 热力耦合符号

热应变进入机械方程时：

```math
\sigma=D(\varepsilon-\varepsilon_{th})
```

所以耦合项是负号：

```math
K_{uT}^e
=
-\int B_u^T D\alpha m N_T\,d\Omega
```

漏掉负号，热膨胀方向会反。

### 11.7 $`T_{ref}`$ 不能忽略

热应变依赖：

```math
T-T_{ref}
```

如果代码直接用 $`T`$，默认就是：

```math
T_{ref}=0
```

若 COMSOL 中 reference temperature 默认不是 0，对比会错。

### 11.8 plane stress / plane strain 选错

薄板通常用 plane stress。

长厚结构截面通常用 plane strain。

如果 MATLAB 和 COMSOL 设置不同，即使代码正确，数值也会不一致。

### 11.9 工程剪应变和张量剪应变混淆

代码里的 $`B`$ 第三行一般对应：

```math
\gamma_{xy}=u_{,y}+v_{,x}
```

所以材料矩阵第三项用：

```math
G=\frac{E}{2(1+\nu)}
```

不要再乘或除 2。

---

## 12. 考前最短记忆版

### 12.1 Tri3

```math
A=\frac12\det
\begin{bmatrix}
1&x_1&y_1\\
1&x_2&y_2\\
1&x_3&y_3
\end{bmatrix}
```

```math
B=
\frac1{2A}
\begin{bmatrix}
b_1&0&b_2&0&b_3&0\\
0&c_1&0&c_2&0&c_3\\
c_1&b_1&c_2&b_2&c_3&b_3
\end{bmatrix}
```

```math
K^e=B^TDB\,tA
```

### 12.2 Q4

```math
N_1=\frac14(1-s)(1-t),\quad
N_2=\frac14(1+s)(1-t)
```

```math
N_3=\frac14(1+s)(1+t),\quad
N_4=\frac14(1-s)(1+t)
```

```math
J=[N_{,s};N_{,t}]\,coords
```

```math
[N_{,x};N_{,y}]=J^{-1}[N_{,s};N_{,t}]
```

```math
K^e=\sum_{g=1}^4B_g^TDB_g\det(J_g)t
```

### 12.3 Thermomechanical

```math
\sigma=D(\varepsilon-\alpha(T-T_{ref})[1,1,0]^T)
```

```math
\begin{bmatrix}
K_{uu}&K_{uT}\\
0&K_{TT}
\end{bmatrix}
\begin{bmatrix}
d\\\theta
\end{bmatrix}
=
\begin{bmatrix}
f_u\\f_T
\end{bmatrix}
```

### 12.4 代码主线

```math
main
\rightarrow
buildBCData
\rightarrow
buildSystem
\rightarrow
assembleSystem
\rightarrow
element
\rightarrow
solveSystem
\rightarrow
postprocess
```

代码框架可以压缩为：

1. main
2. mesh
3. BC data
4. element loop
5. assembly
6. solve
7. postprocess
8. plot
