function [ShapeFunction, DerivativeRDirection, DerivativeSDirection, DerivativeTDirection] = ShapeFunctionAndDerivative(ValueR, ValueS, ValueT)
% This function computes the evaluates value of shape functions for a
% 8-noded element and their derivatives with respect to a given (r, s, t)
% coordinate
% Inputs:
% ValueR: the r coordinate in the parent space
% ValueS: the s coordinate in the parent space
% ValueT: the t coordinate in the parent space
% Outputs:
% ShapeFunction: value of the shape functions
% DerivativeRDirection: derivative of each shape function with respect to
% the r coordinate
% DerivativeSDirection: derivative of each shape function with respect to
% the s coordinate
% DerivativeTDirection: derivative of each shape function with respect to
% the t coordinate

% Compute and store the shape functions
ShapeFunction(1) = 0.125*(1-ValueR)*(1-ValueS)*(1-ValueT);
ShapeFunction(2) = 0.125*(1+ValueR)*(1-ValueS)*(1-ValueT);
ShapeFunction(3) = 0.125*(1+ValueR)*(1+ValueS)*(1-ValueT);
ShapeFunction(4) = 0.125*(1-ValueR)*(1+ValueS)*(1-ValueT);
ShapeFunction(5) = 0.125*(1-ValueR)*(1-ValueS)*(1+ValueT);
ShapeFunction(6) = 0.125*(1+ValueR)*(1-ValueS)*(1+ValueT);
ShapeFunction(7) = 0.125*(1+ValueR)*(1+ValueS)*(1+ValueT);
ShapeFunction(8) = 0.125*(1-ValueR)*(1+ValueS)*(1+ValueT);

% Compute respective derivatives and store them
DerivativeRDirection(1) = -0.125*(1-ValueS)*(1-ValueT);
DerivativeRDirection(2) = 0.125*(1-ValueS)*(1-ValueT);
DerivativeRDirection(3) = 0.125*(1+ValueS)*(1-ValueT);
DerivativeRDirection(4) = -0.125*(1+ValueS)*(1-ValueT);
DerivativeRDirection(5) = -0.125*(1-ValueS)*(1+ValueT);
DerivativeRDirection(6) = 0.125*(1-ValueS)*(1+ValueT);
DerivativeRDirection(7) = 0.125*(1+ValueS)*(1+ValueT);
DerivativeRDirection(8) = -0.125*(1+ValueS)*(1+ValueT);

DerivativeSDirection(1) = -0.125*(1-ValueR)*(1-ValueT);
DerivativeSDirection(2) = -0.125*(1+ValueR)*(1-ValueT);
DerivativeSDirection(3) = 0.125*(1+ValueR)*(1-ValueT);
DerivativeSDirection(4) = 0.125*(1-ValueR)*(1-ValueT);
DerivativeSDirection(5) = -0.125*(1-ValueR)*(1+ValueT);
DerivativeSDirection(6) = -0.125*(1+ValueR)*(1+ValueT);
DerivativeSDirection(7) = 0.125*(1+ValueR)*(1+ValueT);
DerivativeSDirection(8) = 0.125*(1-ValueR)*(1+ValueT);

DerivativeTDirection(1) = -0.125*(1-ValueR)*(1-ValueS);
DerivativeTDirection(2) = -0.125*(1+ValueR)*(1-ValueS);
DerivativeTDirection(3) = -0.125*(1+ValueR)*(1+ValueS);
DerivativeTDirection(4) = -0.125*(1-ValueR)*(1+ValueS);
DerivativeTDirection(5) = 0.125*(1-ValueR)*(1-ValueS);
DerivativeTDirection(6) = 0.125*(1+ValueR)*(1-ValueS);
DerivativeTDirection(7) = 0.125*(1+ValueR)*(1+ValueS);
DerivativeTDirection(8) = 0.125*(1-ValueR)*(1+ValueS);

end