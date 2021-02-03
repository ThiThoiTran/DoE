% Sebastien validation of DOEs part: compute intergration error over a
% given benchmark functions set
%% Name of benchmark functions:

%close all
format long

%% Check MDS errors
% NB_max = 12;
% NB_min = 2;
% NB_tune = 2;
% MDS_error = zeros((NB_max - NB_min) + 1,1);
% MDS_stress = zeros((NB_max - NB_min) + 1,1);
% for i = 7 %1:(NB_max-1)
%     NB = i+1;
%     Neck = countSol(NB);
%     [n_neck, n_dim] = size(Neck);
%     D = zeros(n_neck,n_neck);
%     for k = 1:n_neck
%         for h = k+1:n_neck
%             D(k,h) = sqrt(exp(-NB_tune*(d_neck(Neck(k,:),Neck(k,:))))+exp(-NB_tune*(d_neck(Neck(h,:),Neck(h,:))))-2*exp(-NB_tune*(d_neck(Neck(k,:),Neck(h,:))))); %d_neck(Neck(k,:),Neck(h,:));
%         end
%     end
%     Di = D + D' ;
%     [D1,E] = cmdscale(Di);
%     n_mds = size(D1,2);
%     Distance = zeros(n_neck,n_neck);
%     for l = 1:n_neck
%         for m = l+1:n_neck
%             tmp = 0;
%             for k = 1:n_mds
%                 tmp = tmp + (D1(l,k)-D1(m,k))^2;
%             end
%             Distance(l,m) = sqrt(tmp);
%         end
%     end
%     MDS_error(i) = 2*sum(sum(abs(D - Distance)./(Di+eye(n_neck))))/(n_neck^2);
%     % Compute the stress value
%     MDS_stress(i) = sqrt(sum(sum(D-Distance).^2)/sum(sum(D).^2));
% end
% MDS_error = MDS_error';
% figure(1)
% bar(NB_min:1:NB_max,MDS_error)
% xlabel('NB')
% ylabel('MDScale Relative errors ')
% figure(3)
% bar(NB_min:1:NB_max,MDS_stress)
% xlabel('NB')
% ylabel('MDScale stress values ')
% figure(2)
% imagesc(abs(D - Distance)./(Di + eye(n_neck)))
%  figure(3)



%% Display distances and values of criterion
% NB = 6;
% NC = 1;
% Neck = countSol(NB);
% lower = [0*ones(1,NC) zeros(1,NB)];
% upper = [1*ones(1,NC) ones(1,NB)];
% [n_neck, ~] = size(Neck);
% k = 10;
% nDOE = k*(NC+n_neck);
% n_cont = 20;
% n_take = n_neck;
% X = zeros(n_cont*n_take, NC);
% for i = 1:n_take
%     L = lhsdesign(n_cont, NC);
%     extent = upper - lower;
%     for j=1:NC
%         X((i-1)*n_cont+1:i*n_cont,j) = extent(j)*L(:,j) + repmat(lower(j), n_cont, 1);
%     end
% end
%
%
%     D = zeros(n_neck,n_neck);
%     for k = 1:n_neck
%         for h = k+1:n_neck
%             D(k,h) = sqrt(exp(-NB*(d_neck(Neck(k,:),Neck(k,:))))+exp(-NB*(d_neck(Neck(h,:),Neck(h,:))))-2*exp(-NB*(d_neck(Neck(k,:),Neck(h,:))))); %d_neck(Neck(k,:),Neck(h,:));
%         end
%     end
%     Di = D + D' ;
% % if NB <= 6 && NB>=4
% %     mdsdim = NB - 2;
% % elseif NB<=3
% %     mdsdim = 1;
% % else
% %     mdsdim = NB ;
% % end
% opts = statset('MaxIter',1000,'TolFun',1e-6);
% D1 = cmdscale(Di);%,'Options',opts);
% [~,mdsdim] = size(D1);
% Mat = [kron(D1, ones(n_cont,1)) X];
% % save Rcc.mat Mat
% [pick_MDS] = greedy_cont(NC,NB,Mat,nDOE,n_cont,n_neck,0);
% Mat1 = [kron(Neck, ones(n_cont,1)) X];
% [pick_adapt] = greedy_cont(NC,NB,Mat1,nDOE,n_cont,n_neck,1);
% [pick_adaptnoPen] = greedy_cont(NC,NB,Mat1,nDOE,n_cont,n_neck,2);
% % Return to the original necklaces
% cat_MDS = zeros(nDOE,1);
% cat_adapt = zeros(nDOE,1);
% cat_adaptnoPen = zeros(nDOE,1);
% for i = 1:nDOE
%     for  j = 1:n_neck
%         if isequal(pick_MDS(i,1:mdsdim),D1(j,:))
%             cat_MDS(i) = j;
%             break;
%         end
%      end
% end
% for i = 1:nDOE
%     for  j = 1:n_neck
%         if isequal(pick_adapt(i,1:NB),Neck(j,:))
%             cat_adapt(i) = j;
%             break;
%         end
%      end
% end
% for i = 1:nDOE
%     for  j = 1:n_neck
%         if isequal(pick_adaptnoPen(i,1:NB),Neck(j,:))
%             cat_adaptnoPen(i) = j;
%             break;
%         end
%      end
% end
% num_cat_MDS = zeros(n_neck,1);
% num_cat_adapt = zeros(n_neck,1);
% num_cat_adaptnoPen = zeros(n_neck,1);
% for i = 1:n_neck
%     num_cat_MDS(i) = numel(find(cat_MDS == i));
%     num_cat_adapt(i) = numel(find(cat_adapt == i));
%     num_cat_adaptnoPen(i) = numel(find(cat_adaptnoPen == i));
% end
%
% figure(1)
% subplot(1,3,1)
% bar(num_cat_MDS)
% subplot(1,3,2)
% bar(num_cat_adapt)
% subplot(1,3,3)
% bar(num_cat_adaptnoPen)
% figure(3)
%
%
%
% %% Compute the Xi_T error according to the Greedy paper: note that after operating Greedy-MDS
% % We get back the values of binary, and compute the distribution (average
% % distance) is computed based on d_neck and l2 norm
% pick_MDS = [Neck(cat_MDS,:) pick_MDS(:,mdsdim+1:NC+mdsdim)];
% xi_MDS = sqrt(sum(average_distance(NC,NB,Mat1,Mat1,1)-average_distance(NC,NB,pick_MDS,Mat1,1)).^2);
% xi_adapt = sqrt(sum(average_distance(NC,NB,Mat1,Mat1,1)-average_distance(NC,NB,pick_adapt,Mat1,1)).^2);
% xi_adaptnoPen = sqrt(sum(average_distance(NC,NB,Mat1,Mat1,1)-average_distance(NC,NB,pick_adaptnoPen,Mat1,1)).^2);
% fprintf('Error Xi_T for problem with NB = %d\n',NB)
% fprintf('Greedy-MDS \t\t Greedy-adapted \t\t Greedy-adaptnoPen\n')
% disp([xi_MDS xi_adapt xi_adaptnoPen])
%% Plot illustration of DOEs for proposed methods
NB = 5;
NC = 1;
[~,NB_tune] = len_param(NB,[],[]);
Neck = countSol(NB);
x_l = 0;
x_u = 5;
lower = [x_l*ones(1,NC) zeros(1,NB)];
upper = [x_u*ones(1,NC) ones(1,NB)];
[n_neck, n_dim] = size(Neck);
k = 20;
nDOE = k*(NC*n_neck);
n_cont = 30;
n_take = n_neck;

% Compute the median d_neck
D = zeros(n_neck,n_neck);
for i = 1:n_neck
    for j = i+1:n_neck
        D(i,j) = d_neck(Neck(i,:),Neck(j,:));
    end
end
D = D + D';
D(find(D == 0)) = [];
M_Y = median(D(:));

% R-Random
X_Rand = zeros(nDOE, NC); % matrix load the values
L_rand = rand(nDOE, NC);
extent = upper - lower;
for j = 1:NC
    X_Rand(:,j) = extent(j)*L_rand(:,j) + repmat(lower(j), nDOE, 1);
end
Y_Rand = round(rand(nDOE, NB));
Rand_pick = [Y_Rand X_Rand];
% Rounding LHS
X_LHS = zeros(nDOE, NC); % matrix load the values
L_LHS = lhsdesign(nDOE, NC);
extent = upper - lower;
for j = 1:NC
    X_LHS(:,j) = extent(j)*L_LHS(:,j) + repmat(lower(j), nDOE, 1);
end
Y_LHS = round(lhsdesign(nDOE, NB));
LHS_pick = [Y_LHS X_LHS];


% Sobol picking
Sobol_pick = LHS_cat(NC,NB,Neck,lower,upper,nDOE,0);
% LHS picking
A_LHS_pick = LHS_cat(NC,NB,Neck,lower,upper,nDOE,1);
% Greedy-MDS picking
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
    for h = k+1:n_neck
        D(k,h) = sqrt(exp(-NB_tune*(d_neck(Neck(k,:),Neck(k,:))))+exp(-NB_tune*(d_neck(Neck(h,:),Neck(h,:))))-2*exp(-NB_tune*(d_neck(Neck(k,:),Neck(h,:))))); %d_neck(Neck(k,:),Neck(h,:));
    end
end
D = D + D';
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
NC_tune = (NB_tune*M_Y)/(M_X);

% Target sample
Mat = [kron(D_scale, ones(n_cont,1)) X_scale];
% Run Greedy-MDS algorithm
[pick_MDS] = Greedy(NC,NB,Mat,nDOE,0,NC_tune,NB_tune);
% back to values in X
pick_MDS(:,mdsdim+1:end) = pick_MDS(:,mdsdim+1:end).*(max_Mat_x - min_Mat_x)+ min_Mat_x;
% Adapt mixed kernel directly to Greedy algorithm
Mat1 = [kron(Neck, ones(n_cont,1)) X_scale];
[pick_adapt] = Greedy(NC,NB,Mat1,nDOE,1,NC_tune,NB_tune);
% Back to X values
pick_adapt(:,NB+1:end) = pick_adapt(:,NB+1:end).*(max_Mat_x - min_Mat_x)+ min_Mat_x;
% Return to the original necklaces
cat_MDS = zeros(nDOE,1);
cat_adapt = zeros(nDOE,1);
cat_sobol = zeros(nDOE,1);
cat_A_LHS = zeros(nDOE,1);
cat_Rand = zeros(nDOE,1);
cat_LHS = zeros(nDOE,1);
for i = 1:nDOE
    for  j = 1:n_neck
        if myisrotation(Rand_pick(i,1:NB),Neck(j,:))
            cat_Rand(i) = j;
            break;
        end
    end
end
for i = 1:nDOE
    for  j = 1:n_neck
        if myisrotation(LHS_pick(i,1:NB),Neck(j,:))
            cat_LHS(i) = j;
            break;
        end
    end
end
for i = 1:nDOE
    for  j = 1:n_neck
        if isequal(Sobol_pick(i,1:NB),Neck(j,:))
            cat_sobol(i) = j;
            break;
        end
    end
end
for i = 1:nDOE
    for  j = 1:n_neck
        if isequal(A_LHS_pick(i,1:NB),Neck(j,:))
            cat_A_LHS(i) = j;
            break;
        end
    end
end
for i = 1:nDOE
    for  j = 1:n_neck
        if isequal(pick_MDS(i,1:mdsdim),D_scale(j,:)) % sqrt((pick_MDS(i,1:mdsdim)-D2(j,:)).^2) <= eps    % 
            cat_MDS(i) = j;
            break;
        end
    end
end
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
figure(1)
scatterhist(Rand_pick(:,end), cat_Rand)
figure(2)
scatterhist(LHS_pick(:,end), cat_LHS)
figure(3)
scatterhist(Sobol_pick(:,end), cat_sobol)
figure(4)
scatterhist(A_LHS_pick(:,end),cat_A_LHS)
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
 
% num_cat_MDS = zeros(n_neck,1);
% num_cat_adapt = zeros(n_neck,1);
% num_cat_adaptnoPen = zeros(n_neck,1);
% for i = 1:n_neck
%     num_cat_MDS(i) = numel(find(cat_MDS == i));
%     num_cat_adapt(i) = numel(find(cat_adapt == i));
%     num_cat_adaptnoPen(i) = numel(find(cat_adaptnoPen == i));
% end
% 
% % figure(1)
% % subplot(1,3,1)
% % bar(num_cat_MDS)
% % subplot(1,3,2)
% % bar(num_cat_adapt)
% % subplot(1,3,3)
% % bar(num_cat_adaptnoPen)
% % figure(2)





%% Run different types of DOEs in a set of benchmark functions

n_prob = 15;
name_1 = 'CB2';
name_2 = 'rosen';
name_3 = 'pentagon';
name_4 = 'wong2';
name_5 = 'wong3';
name_6 = 'Branin';
name_7 = 'Hartman3';
name_8 = 'Hartman6';
name_9 = 'perm6';
name_10 = 'perm8';
name_11 = 'Sport';
name_12 = 'ehart3';
name_13 = 'eperm8';
name_14 = 'ebra8';
name_15 = 'ebra10';
%MC_ref = dlmread('MCIE_Vegas_01.txt');
rosen_ref = dlmread('Rosen.txt');
pentagon_ref = dlmread('Pentagon.txt');
CB2_ref = dlmread('CB2.txt');
wong2_ref = dlmread('Wong2.txt');
wong3_ref = dlmread('Wong3.txt');
sport_ref = dlmread('Sport.txt');
Branin_ref = dlmread('Branin.txt');
Hartman3_ref = dlmread('Hartman3.txt');
Hartman6_ref = dlmread('Hartman6.txt');
perm6_ref = dlmread('Perm6.txt');
perm8_ref = dlmread('Perm8.txt');
ebra10_ref = dlmread('EBranin.txt');

M_dim = [2 2;4 3;6 2;10 4;20  6;1 3;2 3; 5 3;5 3;7 3;14 3;2 8;7 8;1 8;1 10];
M_bound = repmat([0 2],n_prob,1);
M_level = [3 4 3 6 14 4 4 4 4 4 4 36 36 36 108] ;
%% Create and compute the error for each function
mode = 2;
name = eval(sprintf('name_%d',mode));
MC_ref = rosen_ref;
NC = eval(sprintf('M_dim(%d,1)',mode));
NB = eval(sprintf('M_dim(%d,2)',mode));
dim = NB + NC;
Neck = countSol(NB);
n_repeat = 20;
k =  1:1:10;
nDOE =  [ k*(NC+size(Neck,1))];%[10 100:100:1000]; %
n_test = length(nDOE);
lower_x = eval(sprintf('M_bound(%d,1)',mode));
upper_x = eval(sprintf('M_bound(%d,2)',mode));
lower = [lower_x*ones(1,NC) zeros(1,NB)];
upper = [upper_x*ones(1,NC) ones(1,NB)];


%% Random picking
eval(sprintf('Random_%s =[];',eval(sprintf('name_%d',mode))));
eval(sprintf('Rand_error_%s = zeros(%d,%d);',eval(sprintf('name_%d',mode)),n_repeat,n_test));
for i = 1:n_test
    tmp = 0;
    err_tmp = 0;
    X = zeros(nDOE(i)*n_repeat, NC); % matrix load the values
    Y = zeros(nDOE(i)*n_repeat, NB);
    for k = 1:n_repeat
        L = rand(nDOE(i), NC);
        extent = upper - lower;
        for j=1:NC
            X((k-1)*nDOE(i)+1:k*nDOE(i),j) = extent(j)*L(:,j) + repmat(lower(j), nDOE(i), 1);
        end
        Y((k-1)*nDOE(i)+1:k*nDOE(i),:) = round(rand(nDOE(i), NB));
        [tmp,ind_rand] = mean_DOE(NC,NB,[Y((k-1)*nDOE(i)+1:k*nDOE(i),:) X((k-1)*nDOE(i)+1:k*nDOE(i),:)],Neck,nDOE(i),mode,M_level(mode));
        err_tmp = eval(sprintf('1/%d * sum(abs(MC_ref(1:%d)-tmp)/abs(MC_ref(1:%d)));',M_level(mode),M_level(mode),M_level(mode)));% sum(abs((Rosen_ref(1:4)-tmp))/abs(Rosen_ref(1:4)));%   abs(MC_ref(mode)-tmp)/abs(MC_ref(mode));
        % eval(sprintf('Rand_error_%s(%d) = err_tmp;',eval(sprintf('name_%d',mode)),k)); % for n_repeat >>
        eval(sprintf('Rand_error_%s(%d,%d) = err_tmp;',eval(sprintf('name_%d',mode)),k,i)); % for nDOE varies
    end
    Rand_pick = [Y X];
    % Compute the mean value of function at this DOEs
    eval(sprintf('Random_%s =[Random_%s ; Rand_pick];',eval(sprintf('name_%d',mode)),eval(sprintf('name_%d',mode))));
end
if n_repeat > 1
     eval(sprintf('Rand_Er_%s = mean(Rand_error_%s);',eval(sprintf('name_%d',mode)),eval(sprintf('name_%d',mode))));
else
    eval(sprintf('Rand_Er_%s = Rand_error_%s;',eval(sprintf('name_%d',mode)),eval(sprintf('name_%d',mode))));
end
    
eval(sprintf('save Rand_Er_%s.mat Rand_Er_%s;',eval(sprintf('name_%d',mode)),eval(sprintf('name_%d',mode))));
eval(sprintf('save Random_%s.mat Random_%s;',eval(sprintf('name_%d',mode)),eval(sprintf('name_%d',mode))));
%% LHS picking
eval(sprintf('LHS_%s =[];',eval(sprintf('name_%d',mode))));
eval(sprintf('LHS_error_%s = zeros(%d,%d);',eval(sprintf('name_%d',mode)),n_repeat,n_test));
for i = 1:n_test
    tmp = 0;
    err_tmp = 0;
    %LHS_pick = zeros(nDOE(i)*n_repeat, NB+NC);
    X = zeros(nDOE(i)*n_repeat, NC); % matrix load the values
    Y = zeros(nDOE(i)*n_repeat, NB);
    for k = 1:n_repeat
        L = lhsdesign(nDOE(i), NC);
        extent = upper - lower;
        for j=1:NC
            X((k-1)*nDOE(i)+1:k*nDOE(i),j) = extent(j)*L(:,j) + repmat(lower(j), nDOE(i), 1);
        end
        Y((k-1)*nDOE(i)+1:k*nDOE(i),:) = round(lhsdesign(nDOE(i), NB));
        [tmp, ind_lhs] = mean_DOE(NC,NB,[Y((k-1)*nDOE(i)+1:k*nDOE(i),:) X((k-1)*nDOE(i)+1:k*nDOE(i),:)],Neck,nDOE(i),mode,M_level(mode));
        err_tmp = eval(sprintf('1/%d * sum(abs(MC_ref(1:%d)-tmp)/abs(MC_ref(1:%d)));',M_level(mode),M_level(mode)-1,M_level(mode)-1));%abs(MC_ref(mode)-tmp)/abs(MC_ref(mode));
        % eval(sprintf('LHS_error_%s(%d) = err_tmp;',eval(sprintf('name_%d',mode)),k));
        % eval(sprintf('LHS_error_%s(%d) = err_tmp;',eval(sprintf('name_%d',mode)),i));
        eval(sprintf('LHS_error_%s(%d,%d) = err_tmp;',eval(sprintf('name_%d',mode)),k,i));
    end
    LHS_pick = [Y X];
    eval(sprintf('LHS_%s =[LHS_%s ; LHS_pick];',eval(sprintf('name_%d',mode)),eval(sprintf('name_%d',mode))));
end
if n_repeat > 1
    eval(sprintf('LHS_Er_%s = mean(LHS_error_%s);',eval(sprintf('name_%d',mode)),eval(sprintf('name_%d',mode))));
else
    eval(sprintf('LHS_Er_%s = LHS_error_%s;',eval(sprintf('name_%d',mode)),eval(sprintf('name_%d',mode))));
end
    
eval(sprintf('save LHS_Er_%s.mat LHS_Er_%s;',eval(sprintf('name_%d',mode)),eval(sprintf('name_%d',mode))));
eval(sprintf('save LHS_%s.mat LHS_%s;',eval(sprintf('name_%d',mode)),eval(sprintf('name_%d',mode))));
%% LHS with the technique using for Sobol (Adapted-LHS)
eval(sprintf('A_LHS_%s =[];',eval(sprintf('name_%d',mode))));
eval(sprintf('A_LHS_error_%s = zeros(%d,%d);',eval(sprintf('name_%d',mode)),n_repeat,n_test));
for i = 1:n_test
    tmp = 0;
    err_tmp = 0;
    A_LHS_pick = zeros(nDOE(i)*n_repeat, NB+NC);
    for k = 1:n_repeat
        A_LHS_pick((k-1)*nDOE(i)+1:k*nDOE(i),:) = LHS_cat(NC,NB,Neck,lower,upper,nDOE(i),'.',1);
       [tmp,ind_A_LHS] = mean_DOE(NC,NB,A_LHS_pick((k-1)*nDOE(i)+1:k*nDOE(i),:),Neck,nDOE(i),mode,M_level(mode));
        err_tmp = eval(sprintf('1/%d * sum(abs(MC_ref(1:%d)-tmp)/abs(MC_ref(1:%d)));',M_level(mode),M_level(mode)-1,M_level(mode)-1));%abs(MC_ref(mode)-tmp)/abs(MC_ref(mode));
       % eval(sprintf('Sobol_error_%s(%d) = err_tmp;',eval(sprintf('name_%d',mode)),k));
        % eval(sprintf('Sobol_error_%s(%d) = err_tmp;',eval(sprintf('name_%d',mode)),i));
        eval(sprintf('A_LHS_error_%s(%d,%d) = err_tmp;',eval(sprintf('name_%d',mode)),k,i));
    end
    eval(sprintf('A_LHS_%s =[A_LHS_%s ; A_LHS_pick];',eval(sprintf('name_%d',mode)),eval(sprintf('name_%d',mode))));
end
if n_repeat > 1
    eval(sprintf('A_LHS_Er_%s = mean(A_LHS_error_%s);',eval(sprintf('name_%d',mode)),eval(sprintf('name_%d',mode))));
else
    eval(sprintf('A_LHS_Er_%s = A_LHS_error_%s;',eval(sprintf('name_%d',mode)),eval(sprintf('name_%d',mode))));
end
eval(sprintf('save A_LHS_Er_%s.mat A_LHS_Er_%s;',eval(sprintf('name_%d',mode)),eval(sprintf('name_%d',mode))));
eval(sprintf('save A_LHS_%s.mat A_LHS_%s;',eval(sprintf('name_%d',mode)),eval(sprintf('name_%d',mode))));


%% Sobol picking
eval(sprintf('Sobol_%s =[];',eval(sprintf('name_%d',mode))));
eval(sprintf('Sobol_error_%s = zeros(%d,%d);',eval(sprintf('name_%d',mode)),n_repeat,n_test));
for i = 1:n_test
    tmp = 0;
    err_tmp = 0;
    Sobol_pick = zeros(nDOE(i)*n_repeat, NB+NC);
    for k = 1:n_repeat
        Sobol_pick((k-1)*nDOE(i)+1:k*nDOE(i),:) = LHS_cat(NC,NB,Neck,lower,upper,nDOE(i),'.',0);
        [tmp,ind_sobol] = mean_DOE(NC,NB,Sobol_pick((k-1)*nDOE(i)+1:k*nDOE(i),:),Neck,nDOE(i),mode,M_level(mode));
        err_tmp = eval(sprintf('1/%d * sum(abs(MC_ref(1:%d)-tmp)/abs(MC_ref(1:%d)));',M_level(mode),M_level(mode)-1,M_level(mode)-1));%abs(MC_ref(mode)-tmp)/abs(MC_ref(mode));
        % eval(sprintf('Sobol_error_%s(%d) = err_tmp;',eval(sprintf('name_%d',mode)),k));
        % eval(sprintf('Sobol_error_%s(%d) = err_tmp;',eval(sprintf('name_%d',mode)),i));
        eval(sprintf('Sobol_error_%s(%d,%d) = err_tmp;',eval(sprintf('name_%d',mode)),k,i));
    end
    eval(sprintf('Sobol_%s =[Sobol_%s ; Sobol_pick];',eval(sprintf('name_%d',mode)),eval(sprintf('name_%d',mode))));
end
if n_repeat > 1
    eval(sprintf('Sobol_Er_%s = mean(Sobol_error_%s);',eval(sprintf('name_%d',mode)),eval(sprintf('name_%d',mode))));
else
    eval(sprintf('Sobol_Er_%s = Sobol_error_%s;',eval(sprintf('name_%d',mode)),eval(sprintf('name_%d',mode))));
end
eval(sprintf('save Sobol_Er_%s.mat Sobol_Er_%s;',eval(sprintf('name_%d',mode)),eval(sprintf('name_%d',mode))));
eval(sprintf('save Sobol_%s.mat Sobol_%s;',eval(sprintf('name_%d',mode)),eval(sprintf('name_%d',mode))));
%% Greedy-MDS picking
eval(sprintf('MDS_%s =[];',eval(sprintf('name_%d',mode))));
eval(sprintf('MDS_error_%s = zeros(%d,%d);',eval(sprintf('name_%d',mode)),n_repeat,n_test));
for i = 1:n_test
    tmp = 0;
    err_tmp = 0;
    MDS_pick = zeros(nDOE(i)*n_repeat, NB+NC);
    for k = 1:n_repeat
        [MDS_pick((k-1)*nDOE(i)+1:k*nDOE(i),:)] = greedy_pick(NC,NB,Neck,lower,upper,nDOE(i),0);
        [tmp,ind_mds] = mean_DOE(NC,NB,MDS_pick((k-1)*nDOE(i)+1:k*nDOE(i),:),Neck,nDOE(i),mode,M_level(mode));
        err_tmp = eval(sprintf('1/%d * sum(abs(MC_ref(1:%d)-tmp)/abs(MC_ref(1:%d)));',M_level(mode),M_level(mode)-1,M_level(mode)-1));%abs(MC_ref(mode)-tmp)/abs(MC_ref(mode));
        % eval(sprintf('MDS_error_%s(%d) = err_tmp;',eval(sprintf('name_%d',mode)),k));
        % eval(sprintf('MDS_error_%s(%d) = err_tmp;',eval(sprintf('name_%d',mode)),i));
        eval(sprintf('MDS_error_%s(%d,%d) = err_tmp;',eval(sprintf('name_%d',mode)),k,i));
    end
    eval(sprintf('MDS_%s =[MDS_%s ; MDS_pick];',eval(sprintf('name_%d',mode)),eval(sprintf('name_%d',mode))));
end
if n_repeat > 1
     eval(sprintf('MDS_Er_%s = mean(MDS_error_%s);',eval(sprintf('name_%d',mode)),eval(sprintf('name_%d',mode))));
else
     eval(sprintf('MDS_Er_%s = MDS_error_%s;',eval(sprintf('name_%d',mode)),eval(sprintf('name_%d',mode))));
end
eval(sprintf('save MDS_Er_%s.mat MDS_Er_%s;',eval(sprintf('name_%d',mode)),eval(sprintf('name_%d',mode))));
eval(sprintf('save MDS_%s.mat MDS_%s;',eval(sprintf('name_%d',mode)),eval(sprintf('name_%d',mode))));
%% Adapted pick
% eval(sprintf('Adapt_%s =[];',eval(sprintf('name_%d',mode))));
% eval(sprintf('Adapt_error_%s = zeros(%d,%d);',eval(sprintf('name_%d',mode)),n_repeat,n_test));
% 
% for i = 1:n_test
%     tmp = 0;
%     err_tmp = 0;
%     Adapt_pick = zeros(nDOE(i)*n_repeat, NB+NC);
%     for k = 1:n_repeat
%         [Adapt_pick((k-1)*nDOE(i)+1:k*nDOE(i),:)] = greedy_pick(NC,NB,Neck,lower,upper,nDOE(i),1);
%         [tmp,ind_adapted] = mean_DOE(NC,NB,Adapt_pick((k-1)*nDOE(i)+1:k*nDOE(i),:),Neck,nDOE(i),mode,M_level(mode));
%         err_tmp = abs(MC_ref(mode)-tmp)/abs(MC_ref(mode));%norm(MC_ref(mode)-tmp,2)^2;
%         % eval(sprintf('Adapt_error_%s(%d) = err_tmp;',eval(sprintf('name_%d',mode)),k));
%         % eval(sprintf('Adapt_error_%s(%d) = err_tmp;',eval(sprintf('name_%d',mode)),i));
%         eval(sprintf('Adapt_error_%s(%d,%d) = err_tmp;',eval(sprintf('name_%d',mode)),k,i));
%     end
%     % Compute the mean value of function at this DOEs
%     eval(sprintf('Adapt_%s =[Adapt_%s ; Adapt_pick];',eval(sprintf('name_%d',mode)),eval(sprintf('name_%d',mode))));
% end
% eval(sprintf('Adapt_Er_%s = mean(Adapt_error_%s);',eval(sprintf('name_%d',mode)),eval(sprintf('name_%d',mode))));
% eval(sprintf('save Adapt_%s.mat Adapt_%s;',eval(sprintf('name_%d',mode)),eval(sprintf('name_%d',mode))));
%% 
% Save Error results
% save Error10repets.mat Rand_error_simulation Rand_Er_simulation LHS_error_simulation LHS_Er_simulation Sobol_error_simulation Sobol_Er_simulation MDS_error_simulation MDS_Er_simulation Adapt_error_simulation Adapt_Er_simulation
fprintf('\n R-Random  \t  R-LHS \t  A-Sobol \t A-LHS \t  MDS \n');
Error = [eval(sprintf('Rand_error_%s(:)',eval(sprintf('name_%d',mode)))) eval(sprintf('LHS_error_%s(:)',eval(sprintf('name_%d',mode))))...
    eval(sprintf('Sobol_error_%s(:)',eval(sprintf('name_%d',mode)))) eval(sprintf('A_LHS_error_%s(:)',eval(sprintf('name_%d',mode))))...
    eval(sprintf('MDS_error_%s(:)',eval(sprintf('name_%d',mode))))];%  eval(sprintf('Adapt_error_%s(:)',eval(sprintf('name_%d',mode))))];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
subplot(1,3,1)
bar(ind_rand(:));
title('Random DOEs')
subplot(1,3,2)
bar(ind_lhs(:));
title('LHS DOEs')
subplot(1,3,3)
bar(ind_A_LHS(:));
title('Adapted-LHS DOEs')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
subplot(1,3,1)
bar(ind_sobol)
title('Sobol sequence DOEs')
subplot(1,3,2)
bar(ind_mds)
title('Greedy-MDS DOEs')
subplot(1,3,3)
bar(ind_A_LHS)
title('A-LHS DOEs')

% figure(3)
% scatter(ones(n_repeat,1),eval(sprintf('Rand_error_%s(:)',eval(sprintf('name_%d',mode)))),'r','filled','jitter','on','jitterAmount',0.1)
% hold on
% scatter(2*ones(n_repeat,1),eval(sprintf('LHS_error_%s(:)',eval(sprintf('name_%d',mode)))),'b','filled','jitter','on','jitterAmount',0.1)
% hold on
% scatter(3*ones(n_repeat,1),eval(sprintf('Sobol_error_%s(:)',eval(sprintf('name_%d',mode)))),'k','filled','jitter','on','jitterAmount',0.1)
% hold on
% scatter(4*ones(n_repeat,1),eval(sprintf('A_LHS_error_%s(:)',eval(sprintf('name_%d',mode)))),'y','filled','jitter','on','jitterAmount',0.1)
% hold on
% scatter(5*ones(n_repeat,1),eval(sprintf('MDS_error_%s(:)',eval(sprintf('name_%d',mode)))),'c','filled','jitter','on','jitterAmount',0.1)
% hold on
% scatter(6*ones(n_repeat,1),eval(sprintf('A_LHS_error_%s(:)',eval(sprintf('name_%d',mode)))),'g','filled','jitter','on','jitterAmount',0.1)
% hold on
% boxplot(Error);
% names = {'Random';  'LHS'; 'Sobol';'Adapted-LHS';'Greedy-MDS';'Greedy-Adapted'};
% set(gca,'xtick',[1:6],'xticklabel',names,'Fontsize',8,'LineWidth',1) %'XTickLabelRotation',45
% ylabel('Integration error')
% %  t = legend('Random','LHS','Sobol','MDS','Adapted')
% %  set(t,'FontSize',8,'LineWidth',1);
%  title('Boxplot and jiter plot for integration error')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
plot(nDOE,eval(sprintf('Rand_Er_%s(:)',eval(sprintf('name_%d',mode)))),'r')
hold on
plot(nDOE,eval(sprintf('LHS_Er_%s(:)',eval(sprintf('name_%d',mode)))),'g')
hold on
plot(nDOE,eval(sprintf('Sobol_Er_%s(:)',eval(sprintf('name_%d',mode)))),'k')
hold on
plot(nDOE,eval(sprintf('MDS_Er_%s(:)',eval(sprintf('name_%d',mode)))),'c')
% hold on
% plot(nDOE,eval(sprintf('Adapt_Er_%s(:)',eval(sprintf('name_%d',mode)))),'r')
hold on
plot(nDOE,eval(sprintf('A_LHS_Er_%s(:)',eval(sprintf('name_%d',mode)))),'y')
legend('Random','LHS','Sobol','Greedy-MDS','Adapted-LHS')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5) % plot in logscale
plot(nDOE,eval(sprintf('log(Rand_Er_%s(:))',eval(sprintf('name_%d',mode)))),'c')
hold on
plot(nDOE,eval(sprintf('log(LHS_Er_%s(:))',eval(sprintf('name_%d',mode)))),'g')
hold on
plot(nDOE,eval(sprintf('log(Sobol_Er_%s(:))',eval(sprintf('name_%d',mode)))),'k')
hold on
plot(nDOE,eval(sprintf('log(MDS_Er_%s(:))',eval(sprintf('name_%d',mode)))),'b')
% hold on
% plot(nDOE,eval(sprintf('log(Adapt_Er_%s(:))',eval(sprintf('name_%d',mode)))),'r')
hold on
plot(nDOE,eval(sprintf('log(A_LHS_Er_%s(:))',eval(sprintf('name_%d',mode)))),'r')
legend('R-Random','R-LHS','A-Sobol','Greedy-MDS','A-LHS')


