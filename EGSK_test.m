% Sebastien validation of DOEs part: compute intergration error over a
% given benchmark functions set
%% illustration
close all
format long




NB = 5;
NC = 1;
Neck = countSol(NB);
x_l = 0;
x_u = 1;
lower = [x_l*ones(1,NC) zeros(1,NB)];
upper = [x_u*ones(1,NC) ones(1,NB)];
[n_neck, n_dim] = size(Neck);
k = 1;
nDOE =  k*(NC*n_neck);
n_cont = 50;
gamma = NB;
if x_l >=0 && x_u <= 1
    scale_x = 0;
else
    scale_x = 1;
end
D = [];
for i = 1:n_neck
    tmp = [];
    for j = i+1:n_neck
        tmp = [tmp d_neck(Neck(i,:),Neck(j,:))];
    end
    D = [D tmp];
end
M_Y = median(D(:));
% Compute the kernel matrix once:
k = zeros(n_neck, n_neck);
for i = 1:n_neck
    for j = i:n_neck
        k(i,j) = string_kernel(Neck(i,:),Neck(j,:),Neck,gamma);
    end
end
K = zeros(n_cont*n_neck, n_cont*n_neck);
for i = 1:n_neck
    for j = i:n_neck
        K((i-1)*n_cont+1:i*n_cont,(j-1)*n_cont+1:j*n_cont) = k(i,j)*ones(n_cont,n_cont);
    end
end
K = triu(K) + triu(K,1)';
%eval(sprintf('save K_NB_%d.mat K;',NB))
%load('K_NB_6.mat')
Distance = zeros(n_neck, n_neck);
for i = 1:n_neck
    for j = i+1:n_neck
        Distance(i,j) = sqrt(string_kernel(Neck(i,:),Neck(i,:),Neck,gamma)+string_kernel(Neck(j,:),Neck(j,:),Neck,gamma)-2*string_kernel(Neck(i,:),Neck(j,:),Neck,gamma));
    end
end
Distance = Distance + (triu(Distance,1))';
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
    NC_tune = len_param(X,gamma,M_Y);
    % Make the big matrix contains a large set of possiblities
   % Mat = [kron(Neck, ones(n_cont,1)) X];
   load('Mat_adapt.mat')
   Mat = Mat(2:end,:);
   [pick] = Greedy_1(NB,Mat,nDOE,Neck,K,gamma,NC_tune, 1);
   F_A_Greedy = subgroup_and_compute_average_f(NC,NB,pick,Neck,4);







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
  %  COmpute the median of X in order to compute the NCtune
   %[NC_tune] = len_param(X,M_Y);
  % Target sample
   % Mat = [kron(D_scale, ones(n_cont,1)) X];
   load('Mat.mat')
   Mat = Mat(2:end,:);
    [pick] = Greedy_1(NB,Mat,nDOE,Neck,K,gamma,[],0);  
    if scale_x
        pick(:,mdsdim+1:end) = pick(:,mdsdim+1:end).*(max_Mat_x - min_Mat_x)+ min_Mat_x;
    end
    % step1: back to the index in MDS
    % step2: from the index, get back the values of y
    cat = zeros(nDOE,1);
    for i = 1:nDOE
        for j = 1:size(D_scale,1)
            if isequal(pick(i,1:mdsdim),D_scale(j,:))
                cat(i) = j;
                break
            end
        end
    end
    pick = [Neck(cat,:) pick(:,mdsdim+1:NC+mdsdim)];





[pick_adapt] = greedy_pick_cluster(NC,NB,Neck,lower,upper,nDOE,M_Y,[],K,gamma,n_cont,scale_x, 1);
cat_adapt = zeros(nDOE,1);
for i = 1:nDOE
    for  j = 1:n_neck
        if isequal(pick_adapt(i,1:NB),Neck(j,:))
            cat_adapt(i) = j;
            break;
        end
    end
end

cat = repmat([1:1:n_neck],n_cont,1);
figure(6)
scatterhist(pick_adapt(:,end),cat_adapt)
figure(7)

%% Check min(eig) for changing NB, gamma
% NB = 2:6;
% j = 1:2:10;
% % Check the minimize eigenvalues
% Min_ei = zeros(length(NB),length(j));
% for i = 1:length(NB)
%     Neck = countSol(NB(i));
%     [n_neck, n_dim] = size(Neck);
%     for j =  1:2:30
%         D = zeros(n_neck,n_neck);
%         for h = 1:n_neck
%             for k = h:n_neck
%                 D(h,k) = string_kernel(Neck(h,:),Neck(k,:),Neck,j) ;
%             end
%         end
%         D = D + (triu(D,1))';
%         Min_ei(i,j) = min(eig(D));
%     end
% end
% %Check min(eig) continuous
% NB = 5:9;
% 
% Min_ei = zeros(length(NB),1);
% for i = 1:length(NB)
%     Neck = countSol(NB(i));
%     [n_neck, n_dim] = size(Neck);
%     lower = [x_l*ones(1,NB(i))];
%     upper = [x_u*ones(1,NB(i))];
%     n_take = 1000;
%     X = zeros(n_take, NB(i));
%     for k = 1:n_take
%         L = lhsdesign(1, NB(i));
%         extent = upper - lower;
%         for j=1:NB(i)
%             X(k,j) = extent(j)*L(:,j) + repmat(lower(j), 1, 1);
%         end
%     end
%     D = zeros(n_neck,n_neck);
%     for h = 1:n_neck
%         for j = h+1:n_neck
%             D(h,j) = exp(-1/2*(norm(X(h,:)-X(j,:),2))^2);%1/2*(norm(X(i,:))+norm(X(j,:)) - norm(X(i,:) - X(j,:)));%
%         end
%     end
%     D = D + D' + eye(n_neck);
%     Min_ei(i) = min(eig(D));  
% end
% disp(Min_ei)

close all
format long
NB = 6;
NC = 1;
Neck = countSol(NB);
x_l = 0;
x_u = 1;
lower = [x_l*ones(1,NB) zeros(1,NB)];
upper = [x_u*ones(1,NB) ones(1,NB)];
[n_neck, n_dim] = size(Neck);
k = 1;
nDOE = k*(NC*n_neck);
n_cont = 30;
n_take = n_neck;
lambda = NB;
 NB_tune = lambda;


% % X = zeros(n_take, NB);
% % for i = 1:n_take
% %     L = lhsdesign(1, NB);
% %     extent = upper - lower;
% %     for j=1:NB
% %         X(i,j) = extent(j)*L(:,j) + repmat(lower(j), 1, 1);
% %     end
% % end

% % D = zeros(n_neck,n_neck);
% % for i = 1:n_neck
% %     for j = i+1:n_neck
% %         D(i,j) = exp(-1/2*(norm(X(i,:)-X(j,:),2))^2);%1/2*(norm(X(i,:))+norm(X(j,:)) - norm(X(i,:) - X(j,:)));%
% %     end
% % end
% % D = D + D' + eye(n_neck);
% % disp(min(eig(D)))



% Compute the median d_neck
% D = zeros(n_neck,n_neck);
% 
% for i = 1:n_neck
%     for j = i+1:n_neck
%         D(i,j) =  d_neck(Neck(i,:),Neck(j,:)); 
%     end
% end
% D = D + D';
% D(find(D == 0)) = [];
% M_Y = 2*median(D(:));

% R-Random
% % X_Rand = zeros(nDOE, NC); % matrix load the values
% % L_rand = rand(nDOE, NC);
% % extent = upper - lower;
% % for j = 1:NC
% %     X_Rand(:,j) = extent(j)*L_rand(:,j) + repmat(lower(j), nDOE, 1);
% % end
% % Y_Rand = round(rand(nDOE, NB));
% % Rand_pick = [Y_Rand X_Rand];
% % % Rounding LHS
% % X_LHS = zeros(nDOE, NC); % matrix load the values
% % L_LHS = lhsdesign(nDOE, NC);
% % extent = upper - lower;
% % for j = 1:NC
% %     X_LHS(:,j) = extent(j)*L_LHS(:,j) + repmat(lower(j), nDOE, 1);
% % end
% % Y_LHS = round(lhsdesign(nDOE, NB));
% % LHS_pick = [Y_LHS X_LHS];
% % 
% % 
% % % Sobol picking
% % Sobol_pick = LHS_cat(NC,NB,Neck,lower,upper,nDOE,0);
% % % LHS picking
% % A_LHS_pick = LHS_cat(NC,NB,Neck,lower,upper,nDOE,1);
% Greedy-MDS picking
%  lambda = 1/M_Y;
% NB_tune = lambda;
X = zeros(n_cont*n_take, NC);
for i = 1:n_take
    L = lhsdesign(n_cont, NC);
    extent = upper - lower;
    for j=1:NC
        X((i-1)*n_cont+1:i*n_cont,j) = extent(j)*L(:,j) + repmat(lower(j), n_cont, 1);
    end
end
D = zeros(n_neck,n_neck);
for k = 1:n_neck
    for h = k:n_neck
        D(k,h) = sqrt(string_kernel(Neck(k,:),Neck(k,:),Neck,lambda)+string_kernel(Neck(h,:),Neck(h,:),Neck,lambda)-2*string_kernel(Neck(k,:),Neck(h,:),Neck,lambda)); %d_neck(Neck(k,:),Neck(h,:));
    end
end
D = D + (triu(D,1))';
[D1,E] = cmdscale(D);
mdsdim = size(D1,2);
fprintf('mds dimension: %d\n\n',mdsdim);
% Scale each dimension to [0 1]
max_Mat = max(D1);
min_Mat = min(D1);
%D_scale = (D1 - min_Mat)./(max_Mat - min_Mat);
D_scale = (D1 - repmat(min_Mat,size(D1,1),1))./(repmat(max_Mat,size(D1,1),1) - repmat(min_Mat,size(D1,1),1));
max_Mat_x = max(X);
min_Mat_x = min(X);
%X_scale = (X - min_Mat_x)./(max_Mat_x - min_Mat_x);
X_scale = (X - repmat(min_Mat_x,size(X,1),1))./(repmat(max_Mat_x,size(X,1),1) - repmat(min_Mat_x,size(X,1),1));


% median x_scale
Mat_X = zeros(size(X_scale,1),size(X_scale,1));
for i = 1:size(X_scale,1)
    for j = i+1:size(X_scale,1)
        Mat_X(i,j) = (norm(X_scale(i,:)-X_scale(j,:)))^2;
    end
end
Mat_X = Mat_X + Mat_X';
Mat_X(find(Mat_X == 0)) = [];
M_X = median(Mat_X(:));
%NC_tune = (NB_tune*M_Y)/(M_X);
NC_tune = (NB_tune)/(M_X*M_Y);


% Target sample
Mat = [kron(D_scale, ones(n_cont,1)) X_scale];
% Run Greedy-MDS algorithm
% [pick_MDS] = Greedy(NC,NB,Mat,nDOE,0,NC_tune,NB_tune);
[pick_MDS] = Greedy_1(NB,Mat,nDOE,Neck,0,lambda,NC_tune, 0);
% back to values in X
pick_MDS(:,mdsdim+1:end) = pick_MDS(:,mdsdim+1:end).*(max_Mat_x - min_Mat_x)+ min_Mat_x;
% Adapt mixed kernel directly to Greedy algorithm
Mat1 = [kron(Neck, ones(n_cont,1)) X_scale];
% [pick_adapt] = Greedy(NC,NB,Mat1,nDOE,1,NC_tune,NB_tune);
[pick_adapt] = greedy_pick_cluster(NC,NB,Neck,lower,upper,nDOE,M_Y,Distance,lambda, 1); 
% Back to X values
pick_adapt(:,NB+1:end) = pick_adapt(:,NB+1:end).*(max_Mat_x - min_Mat_x)+ min_Mat_x;
% Return to the original necklaces
cat_MDS = zeros(nDOE,1);
cat_adapt = zeros(nDOE,1);
cat_sobol = zeros(nDOE,1);
cat_A_LHS = zeros(nDOE,1);
cat_Rand = zeros(nDOE,1);
cat_LHS = zeros(nDOE,1);
% for i = 1:nDOE
%     for  j = 1:n_neck
%         if myisrotation(Rand_pick(i,1:NB),Neck(j,:))
%             cat_Rand(i) = j;
%             break;
%         end
%     end
% end
% for i = 1:nDOE
%     for  j = 1:n_neck
%         if myisrotation(LHS_pick(i,1:NB),Neck(j,:))
%             cat_LHS(i) = j;
%             break;
%         end
%     end
% end
% for i = 1:nDOE
%     for  j = 1:n_neck
%         if isequal(Sobol_pick(i,1:NB),Neck(j,:))
%             cat_sobol(i) = j;
%             break;
%         end
%     end
% end
% for i = 1:nDOE
%     for  j = 1:n_neck
%         if isequal(A_LHS_pick(i,1:NB),Neck(j,:))
%             cat_A_LHS(i) = j;
%             break;
%         end
%     end
% end
% for i = 1:nDOE
%     for  j = 1:n_neck
%         if isequal(pick_MDS(i,1:mdsdim),D_scale(j,:)) % sqrt((pick_MDS(i,1:mdsdim)-D2(j,:)).^2) <= eps    % 
%             cat_MDS(i) = j;
%             break;
%         end
%     end
% end
for i = 1:nDOE
    for  j = 1:n_neck
        if isequal(pick_adapt(i,1:NB),Neck(j,:))
            cat_adapt(i) = j;
            break;
        end
    end
end

cat = repmat([1:1:n_neck],n_cont,1);
%plot illustration
% figure(1)
% scatterhist(Rand_pick(:,end), cat_Rand)
% figure(2)
% scatterhist(LHS_pick(:,end), cat_LHS)
% figure(3)
% scatterhist(Sobol_pick(:,end), cat_sobol)
% figure(4)
% scatterhist(A_LHS_pick(:,end),cat_A_LHS)
figure(5)
scatterhist(pick_MDS(:,end),cat_MDS)
figure(6)
scatterhist(pick_adapt(:,end),cat_adapt)
figure(7)
% scatter(X,cat(:),'b')
% hold on
scatter(pick_MDS(:,end),cat_MDS,'r')
figure(8)
% scatter(X,cat(:),'b')
% hold on
scatter(pick_adapt(:,end),cat_adapt,'r')
figure(9)
 %% Plot stress values 
% % NB_max = 12;
% % NB_min = 2;
% % MDS_error = zeros((NB_max - NB_min) + 1,1);
% % MDS_stress = zeros((NB_max - NB_min) + 1,1);
% % for i = 1:(NB_max-1)
% %     NB = i+1;
% %     lambda = NB;
% %     Neck = countSol(NB);
% %     [n_neck, n_dim] = size(Neck);
% %     D = zeros(n_neck,n_neck);
% %     for k = 1:n_neck
% %         for h = k:n_neck
% %             D(k,h) = sqrt(string_kernel(Neck(k,:),Neck(k,:),Neck,lambda)+string_kernel(Neck(h,:),Neck(h,:),Neck,lambda)-2*string_kernel(Neck(k,:),Neck(h,:),Neck,lambda));
% %         end
% %     end
% %     Di = D + (triu(D,1))' ;
% %     [D1,E] = cmdscale(Di);
% %     n_mds = size(D1,2);
% %     Distance = zeros(n_neck,n_neck);
% %     for l = 1:n_neck
% %         for m = l+1:n_neck
% %             tmp = 0;
% %             for k = 1:n_mds
% %                 tmp = tmp + (D1(l,k)-D1(m,k))^2;
% %             end
% %             Distance(l,m) = sqrt(tmp);
% %         end
% %     end
% %     MDS_error(i) = 2*sum(sum(abs(D - Distance)./(Di+eye(n_neck))))/(n_neck^2);
% %     % Compute the stress value
% %     MDS_stress(i) = sqrt(sum(sum(D-Distance).^2)/sum(sum(D).^2));
% % end
% % MDS_error = MDS_error';
% % figure(1)
% % bar(NB_min:1:NB_max,MDS_error)
% % xlabel('NB')
% % ylabel('MDScale Relative errors ')
% % figure(3)
% % bar(NB_min:1:NB_max,MDS_stress)
% % xlabel('NB')
% % ylabel('MDScale stress values ')
% % figure(4)

% D = zeros(n_neck,n_neck);
% for h = 1:n_neck
%     for k = h:n_neck
%         D(h,k) = string_kernel(Neck(h,:),Neck(k,:),Neck,lambda) ;
%     end
% end
% D = D + (triu(D,1))';
% disp(min(eig(D)));
%% Check min(eig) for changing NB, gamma
% % j = 1:2:30;
% % % Check the minimize eigenvalues
% % Min_ei = zeros(length(NB),length(j));
% % for i = 1:length(NB)
% %     Neck = countSol(NB(i));
% %     [n_neck, n_dim] = size(Neck);
% %     for j = 1:2:30
% %         D = zeros(n_neck,n_neck);
% %         for h = 1:n_neck
% %             for k = h:n_neck
% %                 D(h,k) = string_kernel(Neck(h,:),Neck(k,:),Neck,j) ;
% %             end
% %         end
% %         D = D + (triu(D,1))';
% %         Min_ei(i,j) = min(eig(D));
% %     end
% % end
%% Check min(eig) continuous
% % NB = 5:9;
% % 
% % Min_ei = zeros(length(NB),1);
% % for i = 1:length(NB)
% %     Neck = countSol(NB(i));
% %     [n_neck, n_dim] = size(Neck);
% %     lower = [x_l*ones(1,NB(i))];
% %     upper = [x_u*ones(1,NB(i))];
% %     n_take = 1000;
% %     X = zeros(n_take, NB(i));
% %     for k = 1:n_take
% %         L = lhsdesign(1, NB(i));
% %         extent = upper - lower;
% %         for j=1:NB(i)
% %             X(k,j) = extent(j)*L(:,j) + repmat(lower(j), 1, 1);
% %         end
% %     end
% %     D = zeros(n_neck,n_neck);
% %     for h = 1:n_neck
% %         for j = h+1:n_neck
% %             D(h,j) = exp(-1/2*(norm(X(h,:)-X(j,:),2))^2);%1/2*(norm(X(i,:))+norm(X(j,:)) - norm(X(i,:) - X(j,:)));%
% %         end
% %     end
% %     D = D + D' + eye(n_neck);
% %     Min_ei(i) = min(eig(D));  
% % end
% % disp(Min_ei)