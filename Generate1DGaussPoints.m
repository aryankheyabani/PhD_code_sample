function [Points, Weights] = Generate1DGaussPoints(RequestedNumberofPoints)
% This function is used to generate the gauss points with their
% corresponding weights to be used in a numerical integration in one
% dimension
% Inputs:
% RequestedNumberofPoints: number of Gauss points and weights to be
% generated
% Outputs:
% Points: location of the gauss points
% Weights: weights associated to the gauss points

Points = zeros(RequestedNumberofPoints, 1);
Weights = zeros(RequestedNumberofPoints, 1);

if RequestedNumberofPoints == 1
    Points(1) = 0.0;
    Weights(1) = 2.0;

elseif RequestedNumberofPoints == 2
    Points(1) = -0.577350269189626;
    Points(2) = -Points(1);
    Weights(1) = 1;
    Weights(2) = Weights(1);

elseif RequestedNumberofPoints == 3
    Points(1) = -0.774596669241483;
    Points(2) = 0.0;
    Points(3) = -Points(1);
    Weights(1) = 0.555555555555556;
    Weights(2) = 0.888888888888889;
    Weights(3) = Weights(1);

elseif RequestedNumberofPoints == 4
    Points(1) = -0.861136311594053;
    Points(2) = -0.339981043584856;
    Points(3) = -Points(2);
    Points(4) = -Points(1);
    Weights(1) = 0.347854845137454;
    Weights(2) = 0.652145154862546;
    Weights(3) = Weights(2);
    Weights(4) = Weights(1);

elseif RequestedNumberofPoints == 5
    Points(1) = -0.906179845938664;
    Points(2) = -0.538469310105683;
    Points(3) = 0.0;
    Points(4) = -Points(2);
    Points(5) = -Points(1);
    Weights(1) = 0.236926885056189;
    Weights(2) = 0.478628670499366;
    Weights(3) = 0.568888888888889;
    Weights(4) = Weights(2);
    Weights(5) = Weights(1);

end
end