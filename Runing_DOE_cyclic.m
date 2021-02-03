% Runing implementation for DOE cyclic: using K_soft with lambda = 1, p=1
% Given the input
NB = 6;
NC = 1;
Neck = countSol(NB);
nDOE = 14;
x_l = 0;
x_u = 2;
lambda = 1;
lower = [x_l*ones(1,NC) zeros(1,NB)];
upper = [x_u*ones(1,NC) ones(1,NB)];
[n_neck, n_dim] = size(Neck);

%% Pre-compute the distance matrix

% D = zeros(n_neck,n_neck);
% for i = 1:n_neck
%     for j = i+1:n_neck
%         D(i,j) =  d_neck(Neck(i,:),Neck(j,:));
%     end
% end

D = [];
for i = 1:n_neck
    tmp =[];
    for j = i+1:n_neck
        tmp = [tmp d_neck(Neck(i,:),Neck(j,:))];
    end
    D = [D tmp];
end
M_Y = median(D);
%% Pre-compute the dept between all the points
K = zeros(n_neck,n_neck);
for i = 1:n_neck
    for j = i:n_neck
         K(i,j) = string_kernel(Neck(i,:),Neck(j,:),Neck,1) ;
    end
end

pick_adapt = greedy_pick_cluster(NC,NB,Neck,lower,upper,nDOE,M_Y,[],lambda,K,1);

cat_adapt = zeros(nDOE,1);
for i = 1:nDOE
    for  j = 1:n_neck
        if isequal(pick_adapt(i,1:NB),Neck(j,:))
            cat_adapt(i) = j;
            break;
        end
    end
end
figure(1)
scatterhist(pick_adapt(:,end),cat_adapt)

