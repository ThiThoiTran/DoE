function [pick] = greedy_cont(NC, NB, Rcc,hnum,n_cont,n_neck, mode)
% greedy picking algorithm for picking a required number of points from a
% lagre set of points:
% Input: 
% NB: binary dimension, NC: continuous dimension, Rcc: large matrix
% contains coordinates of points, mode = 0: using the l_2 norm to compute
% the energy criterion, mode = 1: d_mixed = l_2 + d_neck, else dmix=l2+dH
% hnum (= nDOE): number of point requiring (e.g., NB+NC+1 points)
% Return: a matrix contains the set of picking points (hnum*(NB+NC))
% hnum = NB+NC+1;
NC_tune = 50;
NB_tune = 2;
[n_points, n_dim] = size(Rcc);
mdsdim = n_dim - NC;
%% Energycrit
ret = zeros(n_points,1);
if mode == 0% (energycrit) energy is computed according to Eucleadian distance only
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
            currentsum(j) = currentsum(j) + (sqrt(runcrit) );
            
        end
        crit = ret - (currentsum/(i-1));
        [~, id(i)] = min((crit));
        pick(i,:) = Rcc(id(i),:);
    end
elseif mode == 1 % mode = 1: mixed distance: d = d_cont + d_neck
%% energycrit    
    for i = 1:n_points
        runcrit = 0;
        for j = 1:n_points
            runcrittmp = 0;
            for k = NB+1:n_dim
                runcrittmp = runcrittmp + (abs(Rcc(i,k) - Rcc(j,k)))^2;
            end
            runcrit = runcrit + (exp(-NC_tune*runcrittmp)*exp(-sqrt(exp(-NB_tune*(d_neck(Rcc(i,1:NB),Rcc(i,1:NB))))+exp(-NB_tune*(d_neck(Rcc(j,1:NB),Rcc(j,1:NB))))-2*exp(-NB_tune*(d_neck(Rcc(i,1:NB),Rcc(j,1:NB))))))); %(exp(-NB_tune*d_neck(Rcc(i,1:NB), Rcc(j,1:NB)))));
%              runcrit = runcrit + (1/2*(sqrt(norm(Rcc(i,NB+1:n_dim),2))+sqrt(norm(Rcc(i,NB+1:n_dim),2))-sqrt(runcrittmp))*(exp(-1.5*d_neck(Rcc(i,1:NB), Rcc(j,1:NB)))));
       %       runcrit = runcrit + (sqrt(runcrittmp)+sqrt(exp(-NB*(d_neck(Rcc(i,1:NB),Rcc(i,1:NB))))+exp(-NB*(d_neck(Rcc(j,1:NB),Rcc(j,1:NB))))-2*exp(-NB*(d_neck(Rcc(i,1:NB),Rcc(j,1:NB))))));
        end
        ret(i) = runcrit/n_points;
    end
%% Herdingeff    
    id = zeros(hnum,1);
    currentsum = zeros(n_points,1);
    pick = zeros(hnum, n_dim);
    [~,id(1)] = min(ret);
    pick(1,:) = Rcc(id(1),:);
    neck_pen = zeros(n_cont*n_neck,1);
    for i = 2:hnum
        for j = 1:n_points
            runcrit = 0;
            for k = NB+1:n_dim
                runcrit = runcrit + (abs(Rcc(j,k) - pick(i-1,k)))^2;
            end
            currentsum(j) = currentsum(j) + (exp(-NC_tune*runcrit)*exp(-(sqrt(exp(-NB_tune*(d_neck(Rcc(j,1:NB),Rcc(j,1:NB))))+exp(-NB_tune*(d_neck(pick(i-1,1:NB),pick(i-1,1:NB))))-2*exp(-NB_tune*(d_neck(Rcc(j,1:NB),pick(i-1,1:NB))))))));%(exp(-NB_tune*d_neck(Rcc(j,1:NB), pick(i-1,1:NB)))));
%              currentsum(j) = currentsum(j) + (1/2*(sqrt(norm(Rcc(i,NB+1:n_dim),2))+sqrt(norm(pick(i-1,NB+1:n_dim),2))-sqrt(runcrit))*(exp(-1.5*d_neck(Rcc(j,1:NB), pick(i-1,1:NB)))));
 %             currentsum(j) = currentsum(j) + (sqrt(runcrit)+sqrt(exp(-NB*(d_neck(Rcc(j,1:NB),Rcc(j,1:NB))))+exp(-NB*(d_neck(pick(i-1,1:NB),pick(i-1,1:NB))))-2*exp(-NB*(d_neck(Rcc(j,1:NB),pick(i-1,1:NB))))));
        end
        crit = 2*ret - (currentsum/(i-1));
        [~, id(i)] = max((crit));
        pick(i,:) = Rcc(id(i),:);
    end
else % For mixed dneck (no pen) (commented the part of hamming distance and L2 norm)
     for i = 1:n_points
        runcrit = 0;
        for j = 1:n_points
            runcrittmp = 0;
            for k = NB+1:n_dim
                runcrittmp = runcrittmp + (abs(Rcc(i,k) - Rcc(j,k)))^2;
            end
            runcrit = runcrit + (sqrt(runcrittmp)+d_neck(Rcc(i,1:NB), Rcc(j,1:NB)));
        end
        ret(i) = runcrit/n_points;
    end
%% Herdingeff    
    id = zeros(hnum,1);
    val = id;
    currentsum = zeros(n_points,1);
    pick = zeros(hnum, n_dim);
    [~,id(1)] = min(ret);
    pick(1,:) = Rcc(id(1),:);
    neck_pen = zeros(n_cont*n_neck,1);
    for i = 2:hnum
        for j = 1:n_points
            runcrit = 0;
            for k = NB+1:n_dim
                runcrit = runcrit + (abs(Rcc(j,k) - pick(i-1,k)))^2;
            end
            currentsum(j) = currentsum(j) + (sqrt(runcrit)+ d_neck(Rcc(j,1:NB), pick(i-1,1:NB)));
        end
        crit = ret - (currentsum/(i-1));
        [~, id(i)] = min((crit));
        pick(i,:) = Rcc(id(i),:);
    end
    %% Energycrit
%     for i = 1:n_points
%         runcrit = 0;
%         for j = 1:n_points
%             runcrittmp = 0;
%             for k = NB+1:n_dim
%                 runcrittmp = runcrittmp + (abs(Rcc(i,k) - Rcc(j,k)))^2;
%             end
%             runcrit = runcrit + sqrt(runcrittmp)+sum(abs((Rcc(i,1:NB)- Rcc(j,1:NB))));
%         end
%         ret(i) = runcrit/n_points;
%     end
%     %% Herdingeff
%     id = zeros(hnum,1);
%     currentsum = zeros(n_points,1);
%     pick = zeros(hnum, n_dim);
%     [~,id(1)] = min(ret);
%     pick(1,:) = Rcc(id(1),:);
%     for i = 2:hnum
%         crit = zeros(n_points, 1);
%         for j = 1:n_points
%             runcrit = 0;
%             for k = NB+1:n_dim
%                 runcrit = runcrit + (abs(Rcc(j,k) - pick(i-1,k)))^2;
%             end
%             currentsum(j) = currentsum(j) + sqrt(runcrit)+sum(abs((Rcc(i,1:NB)- Rcc(j,1:NB)))) ;
%         end
%         crit = ret - (currentsum/(i-1));
%         [val(i), id(i)] = min(crit);
%         pick(i,:) = Rcc(id(i),:);
%     end
 end
end

