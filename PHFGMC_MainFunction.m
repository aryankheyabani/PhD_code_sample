function [HomogenizedStiffnessMatrix, HomogenizedComplianceMatrix] = PHFGMC_MainFunction(NodalCoordinates, ElementConnectivities, ...
    NumberofElements, FiberProperties, MatrixProperties, FiberorMatrix)

% This is the main function that applies the PHFGMC to compute the
% homogenized material properties based on the constituent material
% properties.
% The function needs to be called as the following:
% [HomogenizedStiffnessMatrix, HomogenizedComplianceMatrix] = PHFGMC_MainFunction(NodalCoordinates, ...
%     ElementConnectivities, NumberofElements, FiberProperties, MatrixProperties, FiberorMatrix);
% Inputs:
% NodalCoordinates: array defining the nodal coordinates of the RUC
% ElementConnectivities: array defining the element connectivities of the
% RUC mesh
% NumberofElements: defines the number of subcells in the RUC
% FiberProperties: array defining the fiber properties
% MatrixProperties: array defining the matrix properties
% FiberorMatrix: array defining whether a subcess is filled with fiber or
% matrix. 1: Fiber, 2: Matrix
% Outputs:
% HomogenizedStiffnessMatrix: stiffness matrix of the composite at the
% macro level
% HomogenizedComplianceMatrix: compliance matrix of the composite at the
% macro level

%% Material Properties

% Get the material properties and assign them to appropriate variables
E11_Fiber = FiberProperties(1);
E22_Fiber = FiberProperties(2);
E33_Fiber = FiberProperties(3);
G12_Fiber = FiberProperties(4);
G13_Fiber = FiberProperties(5);
G23_Fiber = FiberProperties(6);
v12_Fiber = FiberProperties(7);
v13_Fiber = FiberProperties(8);
v23_Fiber = FiberProperties(9);

E_Matrix = MatrixProperties(1);
v_Matrix = MatrixProperties(2);

% Compute stiffness matrices for the fiber and matrix materials
MatrixStiffness = IsotropicStiffnessMatrix(E_Matrix, v_Matrix);
FiberStiffness  = OrthotropicStiffnessMatrix(E11_Fiber, E22_Fiber, ...
    E33_Fiber, G12_Fiber, G23_Fiber, G13_Fiber, v12_Fiber, v23_Fiber,...
    v13_Fiber);

%% Define the facial connectivity arrays
% FacialNodes is an array showing nodes corresponding to each face of an
% element. This array is required to build connectivity between the faces
% of the RUC
% A row of the array is shaped as:
% element number - face number - nodes 1 to 4
FacialNodes = SubcellFacialNodes(NumberofElements, ElementConnectivities);
FacialNodesSize = size(FacialNodes, 1);

% Normally at this point, the FacialNodes array is used in defining the
% connectivity between the subcells of the RUC. Then, an array is formed
% which shows how each face of a subcell is connected with other faces of
% other subcells. This array is called "FacialConnectivities" and is an
% essential part of the code and needs rigorous operations to define it
% based on the nodal and element connectivities of the RUC. However, here
% an already defined array is provided for conciseness of the sample.
% The FacialConnectivies array has the following form in each row:
% subcell a - face a - subcell b - face b
FacialConnectivities = [
    1	1	3	1
    1	3	8	1
    1	4	4	4
    2	1	4	2
    2	3	7	1
    2	4	3	6
    3	2	1	2
    3	3	6	1
    3	4	2	6
    4	1	2	2
    4	3	5	1
    4	6	1	6
    5	2	4	5
    5	3	8	3
    5	4	7	6
    6	2	3	5
    6	3	7	3
    6	4	8	6
    7	2	2	5
    7	5	6	5
    7	4	5	6
    8	2	1	5
    8	5	5	5
    8	4	6	6];

%% Define the general arrays of the system of equations
% GeneralConstants: is the array containing constants required for the
% system of equations
% StrainConstants: contains constants that will be multiplied with the
% remote strains

% Initialize
GeneralConstants = zeros(21*NumberofElements, 21*NumberofElements);
StrainConstants  = zeros(21*NumberofElements, 6);

% Enter displacement continuity and periodicity conditions
dsp1 = 1;
for FaceIndex = 1:size(FacialConnectivities, 1)

    GeneralConstants(3*FaceIndex-2, 21*(FacialConnectivities(FaceIndex,1)-1)+FacialConnectivities(FaceIndex,2)*3-2) = dsp1;
    GeneralConstants(3*FaceIndex-2, 21*(FacialConnectivities(FaceIndex,3)-1)+FacialConnectivities(FaceIndex,4)*3-2) = -dsp1;

    GeneralConstants(3*FaceIndex-1, 21*(FacialConnectivities(FaceIndex,1)-1)+FacialConnectivities(FaceIndex,2)*3-1) = dsp1;
    GeneralConstants(3*FaceIndex-1, 21*(FacialConnectivities(FaceIndex,3)-1)+FacialConnectivities(FaceIndex,4)*3-1) = -dsp1;

    GeneralConstants(3*FaceIndex, 21*(FacialConnectivities(FaceIndex,1)-1)+FacialConnectivities(FaceIndex,2)*3) = dsp1;
    GeneralConstants(3*FaceIndex, 21*(FacialConnectivities(FaceIndex,3)-1)+FacialConnectivities(FaceIndex,4)*3) = -dsp1;

end

% Compute face normals, face areas, and subcell volumes 
% FacialNormals shows the normals to each face of an element:
% element number - face number - unit normal
% FacialAreas shows area of each face for elements:
% element number - face number - area
[FacialAreas, FacialNormals] = FacialNormalsAndAreas(NumberofElements, NodalCoordinates, FacialNodes, FacialNodesSize);

% Compute volume for all the subcells 
NodesperElement = 8;
SubcellVolumes = VolumeCalculation(NumberofElements, NodesperElement,NodalCoordinates, ElementConnectivities);

% Compute and enter traction continuity and periodicity conditions
NumberofEquationsperSubcell = NumberofElements*9;
[GeneralConstants, StrainConstants] = TractionComputations(FacialConnectivities, ElementConnectivities, NodalCoordinates, ...
    FacialNormals, FiberorMatrix, NumberofEquationsperSubcell, FiberStiffness, MatrixStiffness, GeneralConstants, StrainConstants, NodesperElement);

% Compute and store equilibirium related equations
[GeneralConstants, StrainConstants] = EquilibriumComputations(NumberofElements, NodesperElement, ElementConnectivities, NodalCoordinates, FiberorMatrix,...
    FiberStiffness, MatrixStiffness, FacialAreas, FacialNormals, GeneralConstants, StrainConstants, NumberofEquationsperSubcell);

% Compute effective stresses 
[EffectiveStresses, ~, TotalVolume] = EffectiveStressComputations(GeneralConstants, StrainConstants, NodalCoordinates, NodesperElement, FiberorMatrix,...
    FiberStiffness, MatrixStiffness, SubcellVolumes, ElementConnectivities, NumberofElements);

% Compute the homogenized stiffness and compliance matrices 
HomogenizedStiffnessMatrix = EffectiveStresses/TotalVolume;
HomogenizedComplianceMatrix = inv(HomogenizedStiffnessMatrix);
end