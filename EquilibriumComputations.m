function [GeneralConstants, StrainConstants] = EquilibriumComputations(NumberofElements, NodesperElement, ElementConnectivities, NodalCoordinates, FiberorMatrix,...
    FiberStiffness, MatrixStiffness, FacialAreas, FacialNormals, GeneralConstants, StrainConstants, NumberofEquationsperSubcell)
% This function computes equilibrium in subcells of a RUC and assigns it to
% the arrays used to define degrees of freedom for the method.
% The function needs to be called as the following:
% [GeneralConstants, StrainConstants] = EquilibriumComputations(NumberofElements, NodesperElement, ElementConnectivities, NodalCoordinates, FiberorMatrix,...
%     FiberStiffness, MatrixStiffness, FacialAreas, FacialNormals, GeneralConstants, StrainConstants, NumberofEquationsperSubcell);
% Inputs:
% NumberofElements: number of subcells in the RUC
% NodesperElement: number of nodes per each element (subcell)
% ElementConnectivities: connectivity of nodes and the subcells
% NodalCoordinates: defines the coordinates for the nodes in the RUC
% FiberorMatrix: array defining whether each subcell contains fiber or
% matrix
% FiberStiffness: stiffness matrix for the fiber constituent
% MatrixStiffness: stiffness matrix for the matrix constituent
% FacialAreas: array containing the facial area of each subcell
% FacialNormals: array containing the normal for each face
% GeneralConstants: array containing the constants for the system of
% equations
% StrainConstants: array containing constants for D matrix in the system of
% equations (coefficients multiplied with the remote strain vector)
% NumberofEquationsperSubcell: variable defining how many equations of
% equilibirium must be computed for each subcell
% Outputs:
% GeneralConstants: array containing the constants for the system of
% equations
% StrainConstants: array containing constants for D matrix in the system of
% equations (coefficients multiplied with the remote strain vector)

% Define the number of gauss points for each direction
GaussPointsX = 1; GaussPointsY = 1; GaussPointsZ = 1;
[GaussPointLocations, GaussPointWeights] = GenerateGaussPoints(GaussPointsX, GaussPointsY, GaussPointsZ);

AssemblingConstant = NumberofEquationsperSubcell * 2;

% Loop over the elements (subcells), compute the equilibirum equations and
% assemble them to the global system of equations arrays
for ElementIndex = 1:NumberofElements

    % Get nodal coordinates into appropriate arrays conforming to the
    % function inputs used in this loop 
    ElementNodes = zeros(8, 1);
    Xcoordinates = zeros(8, 1);
    Ycoordinates = zeros(8, 1);
    Zcoordinates = zeros(8, 1);

    for i = 1:NodesperElement
        ElementNodes(i) = ElementConnectivities(ElementIndex, i);
        Xcoordinates(i) = NodalCoordinates(ElementNodes(i), 1);
        Ycoordinates(i) = NodalCoordinates(ElementNodes(i), 2);
        Zcoordinates(i) = NodalCoordinates(ElementNodes(i), 3);
    end

    % Get the stiffness matrix of the element by checking whether it is
    % fiber or matrix using the following function
    ElementStiffnessMatrix = GetCurrentStiffnessMatrix(FiberorMatrix, ElementIndex, FiberStiffness, MatrixStiffness);

    % Perform integrations on different faces as well as its contribution
    % to the equlibrium equations
    AuxiliaryArray1 = zeros(3,21);
    AuxiliaryArray2 = zeros(3,6);

    for FaceIndex = 1:6

        for GaussPointXIndex = 1:GaussPointsX

            GaussPointX  = GaussPointLocations(GaussPointXIndex, 1);
            WeightPointX = GaussPointWeights(GaussPointXIndex, 1);

            for GaussPointYIndex = 1:GaussPointsY

                GaussPointY = GaussPointLocations(GaussPointYIndex, 2);
                WeightPointY = GaussPointWeights(GaussPointYIndex, 2);
                
                if FaceIndex == 1
                    [~, DerivativeRDirection, DerivativeSDirection, DerivativeTDirection] = ShapeFunctionAndDerivative(GaussPointX, GaussPointY, -1);
                elseif FaceIndex == 2
                    [~, DerivativeRDirection, DerivativeSDirection, DerivativeTDirection] = ShapeFunctionAndDerivative(GaussPointX, GaussPointY, 1);
                elseif FaceIndex == 3
                    [~, DerivativeRDirection, DerivativeSDirection, DerivativeTDirection] = ShapeFunctionAndDerivative(GaussPointX, -1, GaussPointY);
                elseif FaceIndex == 4
                    [~, DerivativeRDirection, DerivativeSDirection, DerivativeTDirection] = ShapeFunctionAndDerivative(1, GaussPointX, GaussPointY);
                elseif FaceIndex == 5
                    [~, DerivativeRDirection, DerivativeSDirection, DerivativeTDirection] = ShapeFunctionAndDerivative(GaussPointX, 1, GaussPointY);
                elseif FaceIndex == 6
                    [~, DerivativeRDirection, DerivativeSDirection, DerivativeTDirection] = ShapeFunctionAndDerivative(-1, GaussPointX, GaussPointY);
                end
                
                % Compute the jacobian of the element 
                ComputedJacobian = JacobianComputation(NodesperElement, DerivativeRDirection, DerivativeSDirection, DerivativeTDirection, Xcoordinates, Ycoordinates, Zcoordinates);
                % Inverse the jacobian 
                ComputedJacobianInverse = inv(ComputedJacobian);

                % Compute contribution of the gauss point to the integral
                if FaceIndex == 1

                    FaceArea = GetFaceArea(FacialAreas, ElementIndex, FaceIndex);
                    FaceNormal = GetFaceNormal(FacialNormals, ElementIndex, FaceIndex);
                    AuxiliaryArray1 = AuxiliaryArray1 + 0.25*FaceArea*FaceNormal*...
                        ElementStiffnessMatrix*A_mat(ComputedJacobianInverse, GaussPointX, ...
                        GaussPointY, -1)*WeightPointX*WeightPointY;
                
                elseif FaceIndex == 2
  
                    FaceArea = GetFaceArea(FacialAreas, ElementIndex, FaceIndex);
                    FaceNormal = GetFaceNormal(FacialNormals,ElementIndex,FaceIndex);
                    AuxiliaryArray1 = AuxiliaryArray1 + 0.25*FaceArea*FaceNormal*...
                        ElementStiffnessMatrix*A_mat(ComputedJacobianInverse, GaussPointX, ...
                        GaussPointY, 1)*WeightPointX*WeightPointY;


                elseif FaceIndex == 3

                    FaceArea = GetFaceArea(FacialAreas,ElementIndex,FaceIndex);
                    FaceNormal = GetFaceNormal(FacialNormals,ElementIndex,FaceIndex);
                    AuxiliaryArray1 = AuxiliaryArray1 + 0.25*FaceArea*FaceNormal*...
                        ElementStiffnessMatrix*A_mat(ComputedJacobianInverse, GaussPointX,...
                        -1, GaussPointY)*WeightPointX*WeightPointY;

                elseif FaceIndex == 4

                    FaceArea = GetFaceArea(FacialAreas,ElementIndex,FaceIndex);
                    FaceNormal = GetFaceNormal(FacialNormals,ElementIndex,FaceIndex);
                    AuxiliaryArray1 = AuxiliaryArray1 + 0.25*FaceArea*FaceNormal*...
                        ElementStiffnessMatrix*A_mat(ComputedJacobianInverse, 1, GaussPointX,...
                        GaussPointY)*WeightPointX*WeightPointY;

                elseif FaceIndex == 5

                    FaceArea = GetFaceArea(FacialAreas,ElementIndex,FaceIndex);
                    FaceNormal = GetFaceNormal(FacialNormals,ElementIndex,FaceIndex);
                    AuxiliaryArray1 = AuxiliaryArray1 + 0.25*FaceArea*FaceNormal*...
                        ElementStiffnessMatrix*A_mat(ComputedJacobianInverse, GaussPointX,...
                        1, GaussPointY)*WeightPointX*WeightPointY;

                elseif FaceIndex == 6

                    FaceArea = GetFaceArea(FacialAreas,ElementIndex,FaceIndex);
                    FaceNormal = GetFaceNormal(FacialNormals,ElementIndex,FaceIndex);
                    AuxiliaryArray1 = AuxiliaryArray1 + 0.25*FaceArea*FaceNormal*...
                        ElementStiffnessMatrix*A_mat(ComputedJacobianInverse, -1, GaussPointX,...
                        GaussPointY)*WeightPointX*WeightPointY;

                end

            end
        end

        AuxiliaryArray2 = AuxiliaryArray2 + FaceArea*FaceNormal*ElementStiffnessMatrix;

    end
    % Assemble to the general matrix of constants
    GeneralConstants((ElementIndex-1)*3+AssemblingConstant+1:(ElementIndex-1)*3+AssemblingConstant+3, (ElementIndex-1)*21+1:(ElementIndex-1)*21+21) = AuxiliaryArray1;
    StrainConstants((ElementIndex-1)*3+AssemblingConstant+1:(ElementIndex-1)*3+AssemblingConstant+3, :) = -AuxiliaryArray2;

end
end