function ElementStiffnessMatrix = GetCurrentStiffnessMatrix(FiberorMatrix ,ElementIndex, FiberStiffness, MatrixStiffness)
% This function returns stiffness matrix for the current element (subcell)
% based on what material is contained within
% Inputs:
% FiberorMatrix: array defining which elements contain what material
% ElementIndex: the id number of element
% FiberStiffness: stiffness matrix of the fiber
% MatrixStiffness: stiffness matrix of the matrix material
% Outputs:
% ElementStiffnessMatrix: current stiffness matrix for the current element

for Index1 = 1:size(FiberorMatrix, 1)
    if FiberorMatrix(Index1, 1) == ElementIndex
        if FiberorMatrix(Index1, 2) == 1
            ElementStiffnessMatrix = FiberStiffness;
        else
            ElementStiffnessMatrix = MatrixStiffness;
        end
    end
end

end