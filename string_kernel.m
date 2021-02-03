    % COmpute the EGSK:
function [K] = string_kernel(x,y,Neck,gamma)
[n_neck,~] = size(Neck);
% Compute the softmin (integration -> sum discretely)
softmin = 0;
for i = 1:n_neck
     softmin = softmin + exp(-gamma*(d_neck(x,Neck(i,:)) + d_neck(y,Neck(i,:))));
end
K = softmin;
end


