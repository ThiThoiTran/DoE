% function to compute the Of from the average value of x for each level (for both empirical distr and our DOEs)
% Remember the order of M: Y X
function [f] = OF(NC,NB,M,Neck,mode)
n_neck = size(Neck,1);
nDOE = size(M,1); 
x_average = zeros(n_neck,NC);
% Detect all x for each necklace:
for i = 1:n_neck
    X = [];
    for j = 1:nDOE
        if myisrotation(M(j,1:NB),Neck(i,1:NB))
            X = [ X; M(j,NB+1:NB+NC)];
        end
    end
    if ~isempty(X)
         x_average(i,:) = mean(X,1);
    else
         x_average(i,:) = -10*ones(1,NC);
    end
end

% COmpute the average values of OF
F = zeros(n_neck,1);
for i = 1:n_neck
    F(i) = benchmark(mode,x_average(i,:),Neck(i,:),Neck);
end
f = mean(F);
end



            
