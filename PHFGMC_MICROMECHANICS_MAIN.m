%% PHFGMC Main File
% This is the main file used to run the Parametric High Fidelity
% Generalized Method of Cells (PHFGMC) micromechanics.
%
% The PHFGMC is a micromechanics method to obtain the macroscale homogenized
% properties of a periodic composite material based on the properties of its
% constituent materials at the microscale.
%
% To compute the macroscale material properties, a repeating unit cell
% (RUC) is generated based on the periodic structure of the composite
% material. Then, the required computations are performed on the RUC to get
% the material properties. The RUC demonstrates how the constituent
% materials (fiber and matrix in this case) form the microstructure,
% geometries of those material as well as the fiber volume fraction (amount
% of fiber per total volume of the RUC).
%
% So for the PHFGMC, a geometry model is needed to define the
% microstructure. This model is then discretized into subcells similar to
% elements of a FEM RVE (representative volume element) to include the
% different constituent materials. The PHGMC formulation applies different
% boundary conditions on the RUC to get the system of governing equations.
% These include continuity of displacements and tractions between the
% subcells and equilibrium equations in each subcell. Displacement and
% tractions continuities are applied in an average integral sense between
% the faces of subcells inside the RUC. To this end, the code needs those
% face connectivities as input. Herein, the geometry of the RUC is modeled
% in ABAQUS and the face connectivities are already defined and given as
% input to the code. Furthermore, material properties for the constituent
% materials need to be provided.
%
% ********** The code functions in the following steps **********
% 1. The nodal coordinates and element connectivities for the RUC are read
% from the input file
% 2. Material properties for the fiber and matrix are provided
% 3. The code calls the main function to perform PHFGMC calculations and
% generates the compliance matrix of the composite as an output.
% 4. Components of the compliance matrix are then used to define the
% macroscale homogenized material properties such as Modulus of Elasticy E,
% and Shear Modulus G as well as Poisson Ratios.
% * Steps applied in each function are explained within the corresponding
% function.


% Apply the Clear function to free the system memory and clean the
% Workspace
clear;
clc;
warning('off')

% Read the input file
FileIdentifier = fopen('8_elem_small_vf_0.46.txt');

% Enter the number of nodes and elements in the RUC
NumberofNodes = 27;
NumberofElements = 8;

% Initialize the arrays for nodal coordinates and element connectivities
NodalCoordinates = zeros (NumberofNodes,3);
ElementConnectivities = zeros (NumberofElements,8);

% Read data from the input file and assign it to the corresponding arrays
while ~feof(FileIdentifier)
    % Get a line and divide it based on the divider ","
    MyLine = fgetl(FileIdentifier);
    String = split(MyLine, ",");

    % Check if the line starts with "n" and assign it to the nodal
    % coordinates array
    if String(1) == "n"
        i=int16(str2double(String(2)));
        NodalCoordinates(i,1) = str2double(String(3));
        NodalCoordinates(i,2) = str2double(String(4));
        NodalCoordinates(i,3) = str2double(String(5));

        % Check if the line starts with "e" and assign it to the array of
        % element connectivities
    elseif String (1) == "e"
        i=int16(str2double(String(2)));
        ElementConnectivities(i,1) = int16(str2double(String(3)));
        ElementConnectivities(i,2) = int16(str2double(String(4)));
        ElementConnectivities(i,3) = int16(str2double(String(5)));
        ElementConnectivities(i,4) = int16(str2double(String(6)));
        ElementConnectivities(i,5) = int16(str2double(String(7)));
        ElementConnectivities(i,6) = int16(str2double(String(8)));
        ElementConnectivities(i,7) = int16(str2double(String(9)));
        ElementConnectivities(i,8) = int16(str2double(String(10)));
    end
end

% Provide the fiber material properties and assign it to an array so to
% give as an input to the main function. Herein, a glass fiber is assumed
% which has isotropic properties and for this reason values of E, G, and v
% do not change based on the direction. However, the main function requires
% properties for all direction separately so that it can be applied to
% general fiber materials. For this reason, we will provide the same fiber
% properties to all directions separately even though they have the same
% value
E11_Fiber = 72.5e3;
E22_Fiber = 72.5e3;
E33_Fiber = 72.5e3;
v12_Fiber = 0.23;
v23_Fiber = 0.23;
v13_Fiber = 0.23;
G_Fiber = E11_Fiber/(2*(1+v12_Fiber));
G12_Fiber = G_Fiber;  G13_Fiber = G_Fiber; G23_Fiber = G_Fiber;
FiberProperties = [E11_Fiber, E22_Fiber, E33_Fiber, G12_Fiber, ...
    G13_Fiber, G23_Fiber, v12_Fiber, v13_Fiber, v23_Fiber];

% The matrix is always assumed as an isotropic material and for this reason
% it only requires one value for E, G, and v, always.
E_Matrix = 2.9e3; v_Matrix = 0.35;
MatrixProperties = [E_Matrix, v_Matrix];

% The main function also requires which subcell includes fiber and which
% one includes matrix. For this reason, the following array assigns the
% corresponding material to a corresponding subcell number. 1: Fiber, 2:
% Matrix
FiberorMatrix = [1,1;2,2;3,2;4,2;5,2;6,2;7,2;8,1];

% Here, the main PHFGMC function is called to compute the homogenized
% material properties.
% The functions needs to be called as the following:
% [HomogenizedStiffnessMatrix, HomogenizedComplianceMatrix] =
% PHFGMC_MainFunction(...)
% For this example, we only compute the compliance matrix so the stiffness
% matrix output is turned off.
[~, HomogenizedComplianceMatrix] = PHFGMC_MainFunction(NodalCoordinates, ...
    ElementConnectivities, NumberofElements, FiberProperties, MatrixProperties, FiberorMatrix);

% Homogenized material properties for the composite can be obtained by
% considering the components of the compliance matrix as the following.
% The ";" sign is not used to display the computed values in the Command
% Window
E11_Composite = 1/HomogenizedComplianceMatrix(1,1)
E22_Composite = 1/HomogenizedComplianceMatrix(2,2)
E33_Composite = 1/HomogenizedComplianceMatrix(3,3)
G12_Composite = 1/HomogenizedComplianceMatrix(6,6)
G13_Composite = 1/HomogenizedComplianceMatrix(5,5)
G23_Composite = 1/HomogenizedComplianceMatrix(4,4)
v12_Composite = -HomogenizedComplianceMatrix(1,2)*E11_Composite
v13_Composite = -HomogenizedComplianceMatrix(1,3)*E11_Composite
v23_Composite = -HomogenizedComplianceMatrix(2,3)*E22_Composite