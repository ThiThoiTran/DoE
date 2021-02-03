
function [K] = correlation(NC,NB,x,y,theta, mode)
% theta1: lengthscale for continuous variables x
% theta2: lengthscale for discrete variables y
if length(x) ~= length(y)
    printf('check again the input');
end
if mode == 0 % mixed d_Hamming
    K = exp(theta(1)*(-sum(x(NB+1:NB+NC)-y(NB+1:NB+NC))^2) - theta(2)*sum(abs(x(1:NB)-y(1:NB))));
else
    K = exp(theta(1)*(-sum(x(NB+1:NB+NC)-y(NB+1:NB+NC))^2 - theta(2)*d_neck(x(1:NB),y(1:NB))));
end
end