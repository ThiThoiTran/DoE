function [pick] = Greedy(NC, NB, Rcc,hnum, mode)
% greedy picking algorithm
% Input: 
% NB: binary dimension, NC: continuous dimension, Rcc: large matrix to be
% picked contains coordinates of points, hnum: number of DOE's points,
% n_cont: number of continuous points added to one necklace to build Rcc
% n_neck: number of necklaces
% mode = 0: purely continous DOEs
% mode = 1: mixed DOEs with d_mixed = l_2 + d_neck, 
% else: mixed DOEs with dmix=l2+dH 
% Return: a matrix contains the set of picking points (hnum) which
% approximate Rcc distribution
NC_tune = 10;
NB_tune = 2;
[n_points, n_dim] = size(Rcc);
%% Energycrit
ret = zeros(n_points,1);
if mode == 0                 % continuous problems with  Eucleadian distance (Code from Sebastien given in C++)
    for i = 1:n_points
        runcrit = 0;
        for j = 1:n_points
              runcrittmp = 0;
              for k = 1:n_dim
                  runcrittmp = runcrittmp + (abs(Rcc(i,k) - Rcc(j,k)))^2;
              end
              runcrit = runcrit + (sqrt(runcrittmp));
        end
        ret(i) = runcrit/n_points;
    end
    %% Herdingeff
    id = zeros(hnum,1);
    currentsum = zeros(n_points,1);
    pick = zeros(hnum, n_dim);
    [~,id(1)] = min(ret);
    pick(1,:) = Rcc(id(1),:);
    for i = 2:hnum
        for j = 1:n_points
            runcrit = 0;
            for k = 1:n_dim
                runcrit = runcrit + (abs(Rcc(j,k) - pick(i-1,k)))^2;
            end
            currentsum(j) = currentsum(j) + (sqrt(runcrit));
            
        end
        crit = ret - (currentsum/(i));
        [~, id(i)] = min(crit);
        pick(i,:) = Rcc(id(i),:);
    end
elseif mode == 1                % mixed problems with mixed distance: d = d_cont + d_neck
%% energycrit    
    for i = 1:n_points
        runcrit = 0;
        for j = 1:n_points
            runcrittmp = 0;
            for k = NB+1:n_dim
                runcrittmp = runcrittmp + (abs(Rcc(i,k) - Rcc(j,k)))^2;
            end
           % runcrit = runcrit + (exp(-NC_tune*runcrittmp)*(exp(-NB_tune*d_neck(Rcc(i,1:NB), Rcc(j,1:NB)))));
% Sebastien Kernel: 
             runcrit = runcrit + (1/2*((norm(Rcc(i,NB+1:n_dim),2))+(norm(Rcc(j,NB+1:n_dim),2))-sqrt(runcrittmp))*(exp(-NB_tune*(d_neck(Rcc(i,1:NB), Rcc(j,1:NB))))));
        end
        ret(i) = runcrit/n_points;
    end
%% Herdingeff    
    id = zeros(hnum,1);
    currentsum = zeros(n_points,1);
    pick = zeros(hnum, n_dim);
    [~,id(1)] = max(ret);
    pick(1,:) = Rcc(id(1),:);
    for i = 2:hnum
        for j = 1:n_points
            runcrit = 0;
            for k = NB+1:n_dim
                runcrit = runcrit + (abs(Rcc(j,k) - pick(i-1,k)))^2;
            end
          %  currentsum(j) = currentsum(j) + (exp(-NC_tune*runcrit)*(exp(-NB_tune*d_neck(Rcc(j,1:NB), pick(i-1,1:NB)))));
%Sebastien Kernel: 
           currentsum(j) = currentsum(j) + (1/2*((norm(Rcc(i,NB+1:n_dim),2))+ (norm(pick(i-1,NB+1:n_dim),2))-sqrt(runcrit))*(exp(-NB_tune*(d_neck(Rcc(j,1:NB), pick(i-1,1:NB))))));
        end
        crit = (currentsum/(i))-2* ret ;
        [~, id(i)] = min((crit));
        pick(i,:) = Rcc(id(i),:);
    end
else                              % mixed problems with dmixed = Eucledian + dH (commented the part of hamming distance and L2 norm)
    %% Energycrit
    for i = 1:n_points
        runcrit = 0;
        for j = 1:n_points
            runcrittmp = 0;
            for k = NB+1:n_dim
                runcrittmp = runcrittmp + (abs(Rcc(i,k) - Rcc(j,k)))^2;
            end
            runcrit = runcrit + sqrt(runcrittmp)+sum(abs((Rcc(i,1:NB)- Rcc(j,1:NB))));
        end
        ret(i) = runcrit/n_points;
    end
    %% Herdingeff
    id = zeros(hnum,1);
    currentsum = zeros(n_points,1);
    pick = zeros(hnum, n_dim);
    [~,id(1)] = min(ret);
    pick(1,:) = Rcc(id(1),:);
    for i = 2:hnum
        for j = 1:n_points
            runcrit = 0;
            for k = NB+1:n_dim
                runcrit = runcrit + (abs(Rcc(j,k) - pick(i-1,k)))^2;
            end
            currentsum(j) = currentsum(j) + sqrt(runcrit)+sum(abs((Rcc(i,1:NB)- pick(i-1,1:NB)))) ;
        end
        crit = ret - (currentsum/(i));
        [~, id(i)] = min(crit);
        pick(i,:) = Rcc(id(i),:);
    end
 end
end

