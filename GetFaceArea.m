function FaceArea = GetFaceArea(FacialAreas, ElementIndex, FaceIndex)
% This function gets the face area for a face of an element from the
% general array of FacialAreas 
% Inputs: 
% FacialAreas: general array containing areas for the faces of the elements
% ElementIndex: id number of the requested element
% FaceIndex: id number of the requested face of the element 

for Index1 = 1:size(FacialAreas,1)
    if FacialAreas(Index1,1)== ElementIndex && FacialAreas(Index1,2) == FaceIndex
        FaceArea = FacialAreas(Index1,3);
    end
end
end