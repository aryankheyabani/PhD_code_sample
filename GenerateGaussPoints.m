function [GaussPointLocations, GaussPointWeights] = GenerateGaussPoints(GaussPointsX, GaussPointsY, GaussPointsZ)
% This function is used to generate the gauss points with their
% corresponding weights to be used in a numerical integration in three
% dimensions
% Inputs:
% GaussPointsX: number of points in the x direction
% GaussPointsY: number of points in the y direction
% GaussPointsZ: number of points in the z direction
% Outputs:
% GaussPointLocations: location of the gauss points
% GaussPointWeights: weights associated to the gauss points

if GaussPointsX > GaussPointsY
    if GaussPointsX > GaussPointsZ
        NumberofGaussPointLocations = GaussPointsX;
    else
        NumberofGaussPointLocations = GaussPointsZ;
    end
else
    if GaussPointsY > GaussPointsZ
        NumberofGaussPointLocations = GaussPointsY;
    else
        NumberofGaussPointLocations = GaussPointsZ;
    end
end

GaussPointLocations = zeros(NumberofGaussPointLocations, 3);
GaussPointWeights = zeros(NumberofGaussPointLocations, 3);

[PointsX, WeightsX] = Generate1DGaussPoints(GaussPointsX);
[PointsY, WeightsY] = Generate1DGaussPoints(GaussPointsY);
[PointsZ, WeightsZ] = Generate1DGaussPoints(GaussPointsZ);

for IndexX = 1:GaussPointsX
    GaussPointLocations(IndexX,1) = PointsX(IndexX);
    GaussPointWeights(IndexX,1) = WeightsX(IndexX);
end
for IndexY = 1: GaussPointsY
    GaussPointLocations(IndexY,2) = PointsY(IndexY);
    GaussPointWeights(IndexY,2) = WeightsY(IndexY);
end
for IndexZ = 1:GaussPointsZ
    GaussPointLocations(IndexZ,3) = PointsZ(IndexZ);
    GaussPointWeights(IndexZ,3) = WeightsZ(IndexZ);
end
end