%% function to compute 1/Tsum phi(x_i)
function [ret] = average_distance(NC,NB,S,Z,mode)
% Z: ~empirical distribution (Rcc matrix to pick)
% S: our DOEs
[n_points,n_dim] = size(Z);
[n_DOE,~] = size(S);
ret = zeros(n_points,1);
if mode == 0
    for i = 1:n_points
        runcrit = 0;
        for j = 1:n_DOE
            runcrittmp = 0;
            for k = 1:n_dim
                runcrittmp = runcrittmp + (abs(Z(i,k) - S(j,k)))^2;
            end
            runcrit = runcrit + sqrt(runcrittmp);
        end
        ret(i) = runcrit/n_DOE;
    end
else
    for i = 1:n_points
        runcrit = 0;
        for j = 1:n_DOE
            runcrittmp = 0;
            for k = NB+1:n_dim
                runcrittmp = runcrittmp + (abs(Z(i,k) - S(j,k)))^2;
            end
            runcrit = runcrit + (sqrt(runcrittmp)+d_neck(Z(i,1:NB), S(j,1:NB)));
        end
        ret(i) = runcrit/n_DOE;
    end
end
end
    