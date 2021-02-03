function [pick] = greedy_pick_cluster(NC,NB,Neck,lower,upper,hnum,M_Y,Distance,K,lambda,n_cont,scale_x, mode)
n_neck = size(Neck,1);
%% Compute the matrix of distances (d_neck)
if mode == 0
    %% Convert to continuous coordinate for using greedy algorithm
    MDS = cmdscale(Distance);
    mdsdim = size(MDS,2);
    max_Mat = max(MDS);
    min_Mat = min(MDS);
    D_scale = (MDS - min_Mat)./(max_Mat - min_Mat);
    %% Each configuration Y add 100 continuous possiblities
    X = zeros(n_cont*n_neck, NC);
    for i = 1:n_neck
        L = lhsdesign(n_cont, NC);
        extent = upper - lower;
        for j=1:NC
            X((i-1)*n_cont+1:i*n_cont,j) = extent(j)*L(:,j) + repmat(lower(j), n_cont, 1);
        end
    end
    % scale X if scale_x == 1
    if scale_x
        max_Mat_x = max(X);
        min_Mat_x = min(X);
        X = (X - min_Mat_x)./(max_Mat_x - min_Mat_x);
    end
    % COmpute the median of X in order to compute the NCtune
   % [NC_tune] = len_param(X,M_Y);
   % Target sample
    Mat = [kron(D_scale, ones(n_cont,1)) X];
    [pick] = Greedy_1(NB,Mat,hnum,Neck,K,lambda,[],mode);  
    if scale_x
        pick(:,mdsdim+1:end) = pick(:,mdsdim+1:end).*(max_Mat_x - min_Mat_x)+ min_Mat_x;
    end
    % step1: back to the index in MDS
    % step2: from the index, get back the values of y
    cat = zeros(hnum,1);
    for i = 1:hnum
        for j = 1:size(D_scale,1)
            if isequal(pick(i,1:mdsdim),D_scale(j,:))
                cat(i) = j;
                break
            end
        end
    end
    pick = [Neck(cat,:) pick(:,mdsdim+1:NC+mdsdim)];
else
    %% Each configuration Y add 100 continuous possiblities
    X = zeros(n_cont*n_neck, NC);
    for i = 1:n_neck
        L = lhsdesign(n_cont, NC);
        extent = upper - lower;
        for j=1:NC
            X((i-1)*n_cont+1:i*n_cont,j) = extent(j)*L(:,j) + repmat(lower(j), n_cont, 1);
        end
    end
    if scale_x
        max_Mat_x = max(X);
        min_Mat_x = min(X);
        X = (X - min_Mat_x)./(max_Mat_x - min_Mat_x);
    end
    % COmpute the median of X in order to compute the NCtune
    NC_tune = len_param(X,lambda,M_Y);
    % Make the big matrix contains a large set of possiblities
    Mat = [kron(Neck, ones(n_cont,1)) X];
    [pick] = Greedy_1(NB,Mat,hnum,Neck,K,lambda,NC_tune, mode);
    if scale_x
        pick(:,NB+1:NB+NC) = pick(:,NB+1:NB+NC).*(max_Mat_x - min_Mat_x)+ min_Mat_x;
    end
end
end



