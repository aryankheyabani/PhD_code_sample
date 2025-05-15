function FaceNormal = GetFaceNormal(FacialNormals, ElementIndex, FaceIndex)
% This function computes the normal matrix for a face of an element
% Inputs:
% FacialNormals: general array containing the normal vectors for the faces
% of elements
% ElementIndex: id number of the current elment
% FaceIndex: id number of the current face for which the normal matrix is
% requested
% Outputs:
% FaceNormal: the normal matrix for the corresponding requested face

for Index1 = 1:size(FacialNormals,1)
    if FacialNormals(Index1, 1) == ElementIndex && FacialNormals(Index1,2) == FaceIndex
        FaceNormal = [FacialNormals(Index1,3),0,0,0,FacialNormals(Index1,5), FacialNormals(Index1,4);
            0, FacialNormals(Index1,4),0,FacialNormals(Index1,5),0,FacialNormals(Index1,3);
            0,0,FacialNormals(Index1,5),FacialNormals(Index1,4),FacialNormals(Index1,3),0];
    end
end
end