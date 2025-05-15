function [ComputedJacobian] = JacobianComputation(NodesperElement, DerivativeRDirection, DerivativeSDirection, DerivativeTDirection, Xcoordinates, Ycoordinates, Zcoordinates)
% This function computes the three-dimensional jacobian for an element
% based on the given derivatives and coordinates

ComputedJacobian = zeros(3,3);

for i = 1:NodesperElement
    ComputedJacobian(1,1) = ComputedJacobian(1,1) + DerivativeRDirection(i) *Xcoordinates(i);
    ComputedJacobian(1,2) = ComputedJacobian(1,2) + DerivativeRDirection(i) *Ycoordinates(i);
    ComputedJacobian(1,3) = ComputedJacobian(1,3) + DerivativeRDirection(i) *Zcoordinates(i);
    ComputedJacobian(2,1) = ComputedJacobian(2,1) + DerivativeSDirection(i) *Xcoordinates(i);
    ComputedJacobian(2,2) = ComputedJacobian(2,2) + DerivativeSDirection(i) *Ycoordinates(i);
    ComputedJacobian(2,3) = ComputedJacobian(2,3) + DerivativeSDirection(i) *Zcoordinates(i);
    ComputedJacobian(3,1) = ComputedJacobian(3,1) + DerivativeTDirection(i) *Xcoordinates(i);
    ComputedJacobian(3,2) = ComputedJacobian(3,2) + DerivativeTDirection(i) *Ycoordinates(i);
    ComputedJacobian(3,3) = ComputedJacobian(3,3) + DerivativeTDirection(i) *Zcoordinates(i);
end
end