 # FEM Teaching Framework README
# Some code files are provided in only template form, and you need to complete them by yourself.

## 1. Overview

This code package is a lightweight MATLAB finite element framework for teaching and prototyping 2D linear finite element methods with **mixed meshes**. It supports multiple element sets in one model, such as triangular and quadrilateral elements used together in the same analysis. The framework has been organized so that the same workflow can be used for:

- linear elasticity
- steady heat conduction
- linear thermomechanics
- linear piezoelectricity (not required for this course)

The original codes were built around simple element-specific scripts and set-based mesh definitions. In the current framework, that idea is retained and extended so that different element types can coexist in the same model through element sets.

---

## 2. Main features

### 2.1 Mixed mesh support

A model can contain multiple element sets, for example:

- `Tri3_Mech` + `Quad4_Mech`
- `Tri3_Thermal` + `Quad4_Thermal`
- `Tri3_ThermoMech` + `Quad4_ThermoMech`
- `Tri3_Piezo` + `Quad4_Piezo`

Each set is described by:

- the element name
- the number of nodal DoFs
- the number of nodes per element

### 2.2 Unified model workflow

The same overall workflow is used for all physics:

1. define mesh and material data in `model`
2. build the FE system using `buildSystem`
3. build essential BC data using `buildBCData`
4. solve using `solveSystem`
5. postprocess using `extractNodalFields` and `postprocessModel`
6. visualize using the plotting utilities

### 2.3 Available element families

Current element families include:

- **Mechanical**: `Tri3_Mech`, `Quad4_Mech`
- **Thermal**: `Tri3_Thermal`, `Quad4_Thermal`
- **Thermomechanical**: `Tri3_ThermoMech`, `Quad4_ThermoMech`
- **Piezoelectric**: `Tri3_Piezo`, `Quad4_Piezo`

---

## 3. Code structure

A typical code folder contains the following groups of files.

### 3.1 Element routines

These generate element matrices and element vectors.

- `Tri3_Mech.m`
- `Quad4_Mech.m`
- `Tri3_Thermal.m`
- `Quad4_Thermal.m`
- `Tri3_ThermoMech.m`
- `Quad4_ThermoMech.m`
- `Tri3_Piezo.m`
- `Quad4_Piezo.m`

### 3.2 Assembly and solver

- `getDofInfo.m`
- `getElementData.m`
- `getElementDoFs.m`
- `assembleSystem.m`
- `buildSystem.m`
- `buildBCData.m`
- `solveSystem.m`

### 3.3 Natural boundary-condition routines

- `findEdgeOwner.m`
- `Edge2_Thermal.m`
- `assembleNaturalBC_Thermal.m`
- `Edge2_MechTraction.m`
- `assembleNaturalBC_Mech.m`
- `Edge2_ThermoMech.m`
- `assembleNaturalBC_ThermoMech.m`
- `Edge2_Piezo.m`
- `assembleNaturalBC_Piezo.m`

### 3.4 Postprocessing

- `extractNodalFields.m`
- `postprocessModel.m`
- `averageFieldToNodes.m`
- `post_Tri3_Mech.m`
- `post_Quad4_Mech.m`
- `post_Tri3_Thermal.m`
- `post_Quad4_Thermal.m`
- `post_Tri3_ThermoMech.m`
- `post_Quad4_ThermoMech.m`
- `post_Tri3_Piezo.m`
- `post_Quad4_Piezo.m`

### 3.5 Plotting

- `getElementFaces2D.m`
- `plotMesh.m`
- `plotNodalScalar.m`
- `plotElementScalar.m`
- `plotDeformedMesh.m`
- `plotDeformedMeshContour.m`

### 3.6 Validation examples

Suggested example scripts:

- `main_mixed_thermal.m`
- `main_mixed_mech.m`
- `main_mixed_thermomech.m`
- `main_validate_thermal_fourier.m`
- `main_validate_mech_patch.m`
- `main_validate_beam_transverse_refined.m`
- `main_validate_thermomech_free_expansion.m`
- `main_validate_piezo_patch.m`
- `main_validate_piezo_naturalBC_charge.m`

---

## 4. Model definition

All analyses are built through a `model` structure.

### 4.1 Required fields

The most important fields are:

```matlab
model = struct();
model.fieldNames = {...};
model.elesType = {...};
model.elesInfo = struct();
model.elesSetName = {...};
model.nodesInfo = [...];
model.section = [...];
```

### 4.2 `model.fieldNames`

This defines the nodal unknowns of the problem.

Examples:

```matlab
model.fieldNames = {'T'};
model.fieldNames = {'ux','uy'};
model.fieldNames = {'ux','uy','T'};
model.fieldNames = {'ux','uy','phi'};
```

This field is now the main source of truth for DoF interpretation. It is required by the generalized field extraction and system-building functions.

### 4.3 `model.elesType`

This is a cell array defining each element set.

Format of each row:

```matlab
{setID, eleName, nNodeDoF, nEleNodes}
```

Example:

```matlab
model.elesType = {
    1, 'Tri3_Mech',  2, 3;
    2, 'Quad4_Mech', 2, 4
};
```


### 4.4 `model.elesInfo`

This stores the connectivity of each element set.

Example:

```matlab
model.elesInfo.elesSet1 = [
    1, 1, 2, 5;
    2, 1, 5, 4
];

model.elesInfo.elesSet2 = [
    3, 2, 3, 6, 5
];
```

Each row is:

```matlab
[elementID, node1, node2, ...]
```

### 4.5 `model.elesSetName`

This must match the names of the set fields in `model.elesInfo`.

Example:

```matlab
model.elesSetName = {'elesSet1','elesSet2'};
```

It is recommended to define the set names explicitly rather than relying on `fieldnames` order when some sets may be empty.

### 4.6 `model.nodesInfo`

This stores nodal coordinates.

Format:

```matlab
[nodeID, x, y]
```

Example:

```matlab
model.nodesInfo = [
    1, 0.0, 0.0;
    2, 1.0, 0.0;
    3, 2.0, 0.0;
    4, 0.0, 1.0
];
```

### 4.7 `model.section`

This stores the element material/section data. The first column must always be the element ID.

The remaining columns depend on the physics.

Examples:

#### Mechanical

```matlab
[eleID, E, nu]
```

#### Thermal

```matlab
[eleID, kxx, kyy]
```

#### Thermomechanical

```matlab
[eleID, E, nu, alpha, kxx, kyy]
```

#### Piezoelectric

```matlab
[eleID,
 C11 C12 C13 C22 C23 C33,
 e11 e12 e13 e21 e22 e23,
 k11 k12 k22]
```

### 4.8 `model.problem`

Global problem-level data such as thickness, plane stress/strain flag, or reference temperature should be stored here.

Example:

```matlab
model.problem = struct();
model.problem.thickness = 1.0;
model.problem.mechType = 'planeStress';
model.problem.Tref = 0.0;
```

### 4.9 `model.eleLoadData`

This is a cell array with one entry per element.

Example for thermal:

```matlab
model.eleLoadData{e} = struct('heatSource', 0.0);
```

Example for mechanics:

```matlab
model.eleLoadData{e} = struct('bodyForce', [0;0]);
```

Example for piezoelectricity:

```matlab
model.eleLoadData{e} = struct('bodyForce', [0;0], 'chargeDensity', 0.0);
```

### 4.10 `model.nodalVector`

This is the global externally applied nodal vector.

Example:

```matlab
model.nodalVector = zeros(numNodes * nNodeDoF, 1);
```

---

## 5. Boundary conditions

### 5.1 Essential boundary conditions

Essential BCs are defined by:

```matlab
model.essBC = {
    [node list], [value vector];
    [node list], [value vector];
    ...
};
```

Examples:

#### Thermal

```matlab
model.essBC = {
    [1,4,7], [100];
    [3,6,9], [0]
};
```

#### Mechanical

```matlab
model.essBC = {
    [1,4,7], [0,0]
};
```

#### Thermomechanical

```matlab
model.essBC = {
    [1,4,7], [0,0,100]
};
```

#### Piezoelectricity

```matlab
model.essBC = {
    [1,4,7], [0,0,0]
};
```

Essential BC data are converted into prescribed DoFs using `buildBCData.m`.

### 5.2 Natural boundary conditions

Natural BCs are defined by chains of boundary nodes. Adjacent nodes in the chain define the boundary edges. This follows the style already used in the earlier mixed thermal script. fileciteturn1file0

Format:

```matlab
model.natBC = {
    [node1,node2,...], [values];
    ...
};
```

Examples:

#### Thermal

```matlab
model.natBC = {
    [4,1,2,3], [M, S]
};
```

#### Mechanical

```matlab
model.natBC = {
    [3,6,9], [tx, ty]
};
```

#### Thermomechanical

```matlab
model.natBC = {
    [3,6,9], [tx, ty, M, S]
};
```

#### Piezoelectricity

```matlab
model.natBC = {
    [3,6,9], [tx, ty, qs]
};
```

---

## 6. Typical workflow

A standard script should follow this pattern:

```matlab
[globalMatrix, globalVector] = buildSystem(model);
bcData = buildBCData(model);
[reaction, solution] = solveSystem(globalMatrix, globalVector, bcData);

result.u = solution;
result.f = reaction;

nodal = extractNodalFields(model, result.u);
results = postprocessModel(model, result.u);
```

For example, nodal scalar plotting may then be done by:

```matlab
plotNodalScalar(model, nodal.T, 'Temperature');
```

or

```matlab
vm = averageFieldToNodes(model, results, 'vonMises');
plotNodalScalar(model, vm, 'Von Mises Stress');
```

---

## 7. Notes on postprocessing

### 7.1 Nodal field extraction

Use:

```matlab
nodal = extractNodalFields(model, result.u);
```

Returned fields depend on `model.fieldNames`.

Examples:

- thermal: `nodal.T`
- mechanics: `nodal.ux`, `nodal.uy`
- thermomechanics: `nodal.ux`, `nodal.uy`, `nodal.T`
- piezoelectricity: `nodal.ux`, `nodal.uy`, `nodal.phi`

### 7.2 Element-level postprocessing

Use:

```matlab
results = postprocessModel(model, result.u);
```

Each entry contains `results(e).value`, which may include:

- `T`
- `gradT`
- `heatFlux`
- `strain`
- `stress`
- `thermalStrain`
- `Efield`
- `Dfield`
- `vonMises`

### 7.3 Averaging element fields to nodes

Use:

```matlab
nodalField = averageFieldToNodes(model, results, fieldName);
```

This works for both scalar and vector/tensor-like quantities.

Examples:

```matlab
Tnodal = averageFieldToNodes(model, results, 'T');
vm = averageFieldToNodes(model, results, 'vonMises');
stressN = averageFieldToNodes(model, results, 'stress');
qN = averageFieldToNodes(model, results, 'heatFlux');
EN = averageFieldToNodes(model, results, 'Efield');
```

For multi-component quantities, columns correspond to components.

---

## 8. Plotting

### 8.1 Plot mesh

```matlab
plotMesh(model, 'Mesh', true, true);
```

### 8.2 Plot nodal scalar field

```matlab
plotNodalScalar(model, nodal.T, 'Temperature');
```

### 8.3 Plot element scalar field

```matlab
plotElementScalar(model, results, 'vonMises', 'Element von Mises');
```

### 8.4 Plot deformed mesh

```matlab
plotDeformedMesh(model, result.u, 100, 'Deformed Mesh');
plotDeformedMeshContour(model, result.u, 100, 'Deformed Mesh with |u|');
```

These utilities work for mixed Tri3/Quad4 meshes because they use element-set loops and the generalized topology helper `getElementFaces2D.m`.

---

## 9. Validation examples

### 9.1 Thermal validation

Suggested benchmarks:

- linear patch test
- Laplace problem with Fourier-series reference
- manufactured-source Poisson problem

### 9.2 Mechanical validation

Suggested benchmarks:

- constant-strain patch test
- cantilever beam under transverse load

The cantilever beam benchmark is especially useful for comparing pure Quad4, mixed Tri3/Quad4, and pure Tri3 meshes. A typical result is that the reaction is captured accurately while bending deflection is underestimated for Tri3-dominated meshes, which is physically consistent with the stiffness behavior of linear constant-strain triangles.

### 9.3 Thermomechanical validation

Suggested benchmarks:

- free thermal expansion
- restrained thermal expansion

The free-expansion benchmark is especially useful because the exact stress is zero and the exact displacement is linear.

### 9.4 Piezoelectric validation

Suggested benchmarks:

- piezoelectric patch test with exact linear `ux`, `uy`, and `phi`
- electrostatic strip with prescribed potential or surface charge

---

## 10. Common errors and troubleshooting

### 10.1 `model.fieldNames must be defined`

This occurs when `extractNodalFields`, `buildSystem`, or other generalized functions are used without defining the field layout.

Fix by adding one of:

```matlab
model.fieldNames = {'T'};
model.fieldNames = {'ux','uy'};
model.fieldNames = {'ux','uy','T'};
model.fieldNames = {'ux','uy','phi'};
```

### 10.2 Passing `'u'` instead of the solution vector

Incorrect:

```matlab
nodal = extractNodalFields(model, 'u');
```

Correct:

```matlab
nodal = extractNodalFields(model, result.u);
```

### 10.3 Empty element sets

If one element set is empty, keep the corresponding `elesSet` field explicitly defined as an empty array and skip it during assembly/postprocessing.

Example:

```matlab
model.elesInfo.elesSet1 = triElems;
model.elesInfo.elesSet2 = quadElems;
model.elesSetName = {'elesSet1','elesSet2'};
```

### 10.4 Wrong mixed mesh connectivity

If the global matrix has entire zero rows/columns, check whether some nodes are not connected to any element. This happened in an earlier mixed mechanical example where a support node was defined but did not belong to the actual mesh, producing zero rows in the global stiffness.

### 10.5 Piezoelectric plotting fails on `Tri3_Piezo` or `Quad4_Piezo`

Make sure `getElementFaces2D.m` recognizes those names, or better, use a prefix-based version:

```matlab
if startsWith(eleName, 'Tri3')
    ...
elseif startsWith(eleName, 'Quad4')
    ...
end
```

---

## 11. Design philosophy and limitations

### 11.1 Teaching-oriented design

This code package is intended for teaching and lightweight research prototyping, not as a production FE code. The priority is clarity, consistency, and extensibility.

### 11.2 One nodal field layout per model

Within a single model, all element sets must share the same nodal DoF layout. For example, a model may mix:

- `Tri3_Mech` and `Quad4_Mech`

but not:

- `Tri3_Mech` and `Quad4_Thermal`

unless the framework is extended to support field-superset formulations.

### 11.3 Low-order element behavior

Patch tests may pass exactly while bending-dominated problems still show element-dependent stiffness effects. For example, linear constant-strain triangles can be noticeably too stiff in cantilever-bending benchmarks, even when assembly and BC handling are correct.

---

## 12. Minimal working example

```matlab
model = struct();
model.fieldNames = {'T'};

model.elesType = {
    1, 'Tri3_Thermal',  1, 3;
    2, 'Quad4_Thermal', 1, 4
};

model.elesInfo.elesSet1 = [1,1,2,5];
model.elesInfo.elesSet2 = [2,2,3,6,5];
model.elesSetName = {'elesSet1','elesSet2'};

model.nodesInfo = [
    1,0,0;
    2,1,0;
    3,2,0;
    4,0,1;
    5,1,1;
    6,2,1
];

model.problem.thickness = 1.0;
model.section = [
    1,10,10;
    2,10,10
];

model.eleLoadData = {
    struct('heatSource',0.0);
    struct('heatSource',0.0)
};

model.nodalVector = zeros(size(model.nodesInfo,1),1);

model.essBC = {
    [1,4], [100];
    [3,6], [0]
};

model.natBC = {};

[globalMatrix, globalVector] = buildSystem(model);
bcData = buildBCData(model);
[reaction, solution] = solveSystem(globalMatrix, globalVector, bcData);

result.u = solution;
result.f = reaction;

nodal = extractNodalFields(model, result.u);
results = postprocessModel(model, result.u);
plotNodalScalar(model, nodal.T, 'Temperature');
```

---
