% Comparing according to the main article of kernel herding
%% Set up values
close all
NB = 12;
NC = 1;
load('Neck.mat')
%Neck = countSol(NB);
x_l = -0.03;
x_u = 0.03;
scale_x = 1;
lower = [x_l*ones(1,NC) zeros(1,NB)];
upper = [x_u*ones(1,NC) ones(1,NB)];
[n_neck, n_dim] = size(Neck);
h_len = 2;
k = 3:h_len:7;
nDOE = k*NC*n_neck;
n_long = length(k); % equal length of k
n_cont = 10;
n_run_ref = 1; % for runing random DOEs
%mode = 10; % for call benchmark

lambda = NB;
%% Complement computation
load('precompute.mat')
% Distance = zeros(n_neck, n_neck);
% for i = 1:n_neck
%     for j = i:n_neck
%         Distance(i,j) = sqrt(string_kernel(Neck(i,:),Neck(i,:),Neck,lambda)+string_kernel(Neck(j,:),Neck(j,:),Neck,lambda)-2*string_kernel(Neck(i,:),Neck(j,:),Neck,lambda));
%     end
% end
% Distance = Distance + (triu(Distance,1))';
% 
% % Compute the median d_neck
%  D = [];
% for i = 1:n_neck
%     tmp = [];
%     for j = i+1:n_neck
%         tmp = [tmp d_neck(Neck(i,:),Neck(j,:))];
%     end
%     D = [D tmp];
% end
% M_Y = median(D(:));
% k = zeros(n_neck, n_neck);
% for i = 1:n_neck
%     for j = i:n_neck
%         k(i,j) = string_kernel(Neck(i,:),Neck(j,:),Neck,lambda);
%     end
% end
% K = zeros(n_cont*n_neck, n_cont*n_neck);
% for i = 1:n_neck
%     for j = i:n_neck
%         K((i-1)*n_cont+1:i*n_cont,(j-1)*n_cont+1:j*n_cont) = k(i,j)*ones(n_cont,n_cont);
%     end
% end
% K = triu(K) + triu(K,1)';
%% Target distribution (need X for each Neck)
% n_tar = 100;
% F_ref = [];
% for i = 1:n_neck
%     X = zeros(n_tar, NC);
%     L = rand(n_tar, NC);
%     extent = upper - lower;
%     for j=1:NC
%         X(:,j) = extent(j)*L(:,j) + repmat(lower(j), n_tar, 1);
%     end
%     F_ref_tmp = zeros(1,n_tar);
%     for k = 1:n_tar
%         WriteData4R(X(k,:), Neck(i,:));
%         [status,lastEGODisplays] = system('"C:\\Program Files\\R\\R-4.0.1\\bin\\Rscript.exe" --vanilla .\\All_combination.R');
%          F_ref_tmp(k) = ReadDataFromR;
%     end
%     F_ref = [F_ref mean(F_ref_tmp)];
%     
% end
% eval(sprintf('save F_ref.mat F_ref;'))
load('F_ref.mat')


Error_LHS_long = zeros(n_run_ref,n_long);
Error_A_LHS_long = zeros(n_run_ref,n_long);
Error_MDS_long = zeros(n_run_ref,n_long);
Error_A_Greedy_long = zeros(n_run_ref,n_long);



%% Run R-Random and R-LHS 50 times for benchmarks
for i = 1:n_run_ref
    F_LHS_tmp = [];
    F_A_LHS_tmp = [];
    F_MDS_tmp = [];
    F_A_Greedy_tmp = [];
    for h = [2:3]
        % Rounding LHS
        X_LHS = zeros(nDOE(h), NC); % matrix load the values
        L = lhsdesign(nDOE(h), NC);
        extent = upper - lower;
        for j = 1:NC
            X_LHS(:,j) = extent(j)*L(:,j) + repmat(lower(j), nDOE(h), 1);
        end
        Y_LHS = round(lhsdesign(nDOE(h), NB));
        Pick_LHS = [Y_LHS X_LHS];
        %% 
        % LHS picking
        Pick_A_LHS = LHS_cat(NC,NB,Neck,lower,upper,nDOE(h),1);
        % Greedy-MDS
        Pick_MDS = greedy_pick_cluster(NC,NB,Neck,lower,upper,nDOE(h),M_Y,Distance,K,lambda,n_cont,scale_x, 0);
       % Adapted-Greedy
        Pick_A_Greedy = greedy_pick_cluster(NC,NB,Neck,lower,upper,nDOE(h),M_Y,Distance,K,lambda,n_cont,scale_x, 1);
        % Compute the subgroup and average of x in each method
        %subgroup_and_compute_average_safran(NC,NB,Mat,Neck)
        F_LHS = subgroup_and_compute_average_safran(NC,NB,Pick_LHS,Neck);
        %%
        F_A_LHS = subgroup_and_compute_average_safran(NC,NB,Pick_A_LHS,Neck);
        F_MDS = subgroup_and_compute_average_safran(NC,NB,Pick_MDS,Neck);
        F_A_Greedy = subgroup_and_compute_average_safran(NC,NB,Pick_A_Greedy,Neck);
        
        eval(sprintf('save F_LHS_safran_%d.mat F_LHS',nDOE(h)));
        eval(sprintf('save F_A_LHS_safran_%d.mat F_A_LHS',nDOE(h)))
        eval(sprintf('save F_Greedy_MDS_safran_%d.mat F_MDS',nDOE(h)))
        eval(sprintf('save F_A_Greedy_safran_%d.mat F_A_Greedy ',nDOE(h)))
 
        F_LHS_tmp = [F_LHS_tmp F_LHS];
        F_A_LHS_tmp = [F_A_LHS_tmp F_A_LHS];
        F_MDS_tmp = [F_MDS_tmp F_MDS];
        F_A_Greedy_tmp = [F_A_Greedy_tmp F_A_Greedy];
       
        % Compute the errors: compute at each levels, then mean (or median errors of these)
        Err_A_LHS =  (norm(F_ref - F_A_LHS))^2;
        Err_LHS = (norm(F_ref - F_LHS))^2;
        Err_MDS = (norm(F_ref - F_MDS))^2;
        Err_A_Greedy =  (norm(F_ref - F_A_Greedy))^2;
        
        
        Error_LHS_long(i,h) =  sqrt((1/n_neck)*Err_LHS);
        %%
        Error_A_LHS_long(i,h) =  sqrt((1/n_neck)*Err_A_LHS);
        Error_MDS_long(i,h) =  sqrt((1/n_neck)*Err_MDS);
        Error_A_Greedy_long(i,h) =  sqrt((1/n_neck)*Err_A_Greedy);
            
    end

        
end



LHS_mean = median(Error_LHS_long);
A_LHS_mean = median(Error_A_LHS_long);
MDS_mean = median(Error_MDS_long);
A_Greedy_mean = median(Error_A_Greedy_long);




eval(sprintf(' save LHS_soft_safran.mat Error_LHS_long;'))
eval(sprintf(' save A_LHS_soft_safran.mat Error_A_LHS_long;'))
eval(sprintf(' save MDS_soft_safran.mat Error_MDS_long;'))
eval(sprintf(' save A_Greedy_soft_safran.mat Error_A_Greedy_long;'))

eval(sprintf(' save F_LHS_soft_safran.mat F_LHS_tmp;'));
eval(sprintf(' save F_A_LHS_soft_safran.mat F_A_LHS_tmp;'));
eval(sprintf(' save F_MDS_soft_safran.mat F_MDS_tmp;'));
eval(sprintf(' save F_A_Greedy_soft_safran.mat F_A_Greedy_tmp;'));
lwd = 2;
figure(1)
semilogy(nDOE(1:h_len:n_long),A_LHS_mean(1:h_len:n_long),'g','LineWidth',lwd)
hold on
semilogy(nDOE(1:h_len:n_long),MDS_mean(1:h_len:n_long),'b','LineWidth',lwd)
hold on
semilogy(nDOE(1:h_len:n_long),A_Greedy_mean(1:h_len:n_long),'r','LineWidth',lwd)
legend('R-LHS','A-LHS','Greedy-MDS','Adapted-Greedy')


function [OF] = ReadDataFromR
fid=fopen('DataFromR.txt','r');
 OF=fscanf(fid,'%f',[1,1]);
%  OF = mean(V);
fclose(fid);
end
 
function WriteData4R( X, Y)
fid=fopen('DataM.txt','w');
fprintf(fid, '%f ',X);fprintf(fid, '\n');

for i=1: size(Y,2)
    fprintf(fid, '%f ',Y(i));
end
fclose(fid);
end
 
 	
