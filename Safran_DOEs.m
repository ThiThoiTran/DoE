% DOEs for Safan's simulation: benchmark: Wong2
n_repet = 100;
NC = 1;% 1;
NB = 12;%12;
freq_min = -0.03;% [0 1 ]for Hartman3
freq_max = 0.03;
x_scale = 1;
lower = [freq_min*ones(NC,1); zeros(NB,1)];
upper = [freq_max*ones(NC,1); ones(NB,1)];
nDOE = 14;
load('Neck.mat')
Neck = Neck(2:end-1,:);
% Neck = countSol(NB);
n_neck = size(Neck,1);

lambda = n_neck;% NB;
n_cont = 30;
D = [];
for i = 1:n_neck
    tmp = [];
    for j = i+1:n_neck
        tmp = [tmp d_neck(Neck(i,:), Neck(j,:))];
    end
    D = [D tmp];
end
M_Y = median(D);
% Compute matrix of distance
Distance = zeros(n_neck, n_neck);
for i = 1:n_neck
    for j = i+1:n_neck
        Distance(i,j) = sqrt(string_kernel(Neck(i,:),Neck(i,:),Neck,lambda)+string_kernel(Neck(j,:),Neck(j,:),Neck,lambda)-2*string_kernel(Neck(i,:),Neck(j,:),Neck,lambda));
    end
end
Distance = Distance + (triu(Distance,1))';

k = zeros(n_neck, n_neck);
for i = 1:n_neck
    for j = i:n_neck
        k(i,j) = string_kernel(Neck(i,:),Neck(j,:),Neck,lambda);
    end
end
K = zeros(n_cont*n_neck, n_cont*n_neck);
for i = 1:n_neck
    for j = i:n_neck
        K((i-1)*n_cont+1:i*n_cont,(j-1)*n_cont+1:j*n_cont) = k(i,j)*ones(n_cont,n_cont);
    end
end
K = triu(K) + triu(K,1)';
mode = 4;
S = zeros(nDOE*n_repet,NB+NC);
%% 
for i = 1:n_repet
    if mode == 1 % Rounding LHS
        X = zeros(nDOE, NC); % matrix load the values
        L = lhsdesign(nDOE, NC);
        extent = upper - lower;
        for j = 1:NC
            X(:,j) = extent(j)*L(:,j) + repmat(lower(j), nDOE, 1);
        end
        Y = round(lhsdesign(nDOE, NB));
        pick = [Y X];
        
    elseif mode == 2 % Adapted LHS
        pick = LHS_cat(NC,NB,Neck,lower,upper,nDOE,1);
    elseif mode == 3 % Adapted greedy
        pick = greedy_pick_cluster(NC,NB,Neck,lower,upper,nDOE,M_Y,[],K,lambda,n_cont,x_scale, 1);
    elseif mode == 4 % Greedy MDS
        pick = greedy_pick_cluster(NC,NB,Neck,lower,upper,nDOE,M_Y,Distance,K,lambda,n_cont,x_scale, 0);
    end
    S((i-1)*nDOE+1:i*nDOE,:)= pick;
end
  %  save temporary.mat S
 save interpolation_Safran_MDS_14_100rep.mat S  

% L = lhsdesign(n_neck, NC);
% extent = n_neck - 1;
% take = extent*L + ones(n_take, 1);
% Take = Neck(round(take),:);
n_cont = 10;
load('precompute.mat')
% D = [];
% for i = 1:n_take
%     tmp = [];
%     for j = i+1:n_take
%         tmp = [tmp d_neck(Take(i,:), Take(j,:))];
%     end
%     D = [D tmp];
% end
% M_Y = median(D);
% k = zeros(n_take, n_take);
% for i = 1:n_take
%     for j = i:n_take
%         k(i,j) = string_kernel(Take(i,:),Take(j,:),Take,lambda);
%     end
% end
% K = zeros(n_cont*n_take, n_cont*n_take);
% for i = 1:n_take
%     for j = i:n_take
%         K((i-1)*n_cont+1:i*n_cont,(j-1)*n_cont+1:j*n_cont) = k(i,j)*ones(n_cont,n_cont);
%     end
% end
% K = triu(K) + triu(K,1)';

S = zeros(nDOE*n_repet,NB+NC);
for h = 1:n_repet
    % pick = greedy_pick_cluster(NC,NB,Neck,lower,upper,nDOE,M_Y,[],K,lambda,n_cont,x_scale, 1);
    % pick = LHS_cat(NC,NB,Neck,lower,upper,nDOE,1);
    X = zeros(nDOE, NC); % matrix load the values
    L = lhsdesign(nDOE, NC);
    extent = upper - lower;
    for j = 1:NC
        X(:,j) = extent(j)*L(:,j) + repmat(lower(j), nDOE, 1);
    end
    Y = round(lhsdesign(nDOE, NB));
    pick = [Y X];
    S((h-1)*nDOE+1:h*nDOE,:) = pick;
end


save interpolation_Safran_R_LHS_50_5rep.mat S


%% Check if there is any duplicated points
check = 1;
t = zeros(nDOE,1);
for i = 1:nDOE
    k = zeros(i-1,1);
    for j = 1:i-1
        if ~myisrotation(S(i,1:NB),S(j,1:NB))
            k(j) = 1;
        end
    end
    if isequal(k,ones(i-1,1))
        t(i) = 1;
    end
end
if sum(t) == nDOE
    check = 0;
end
check


