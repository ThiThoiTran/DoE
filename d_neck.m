% function compute d_neck
function [d] = d_neck(x,y)
n = length(x);
if n ~= length(y)
    fprintf('Check again the size of input');
end
F = zeros(n,1);
for i = 1:n
    F(i) = sum(abs(x - circshift(y,i)));%[0,i]
end
d = min(F);
end