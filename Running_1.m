% Comparing according to the main article of kernel herding
%% Set up values
close all
NB = 7;
NC = 1;
Neck = countSol(NB);
x_l = 0;
x_u = 3;
scale_x = 1;
lower = [x_l*ones(1,NC) zeros(1,NB)];
upper = [x_u*ones(1,NC) ones(1,NB)];
[n_neck, n_dim] = size(Neck);
k = 5;%1:30;
nDOE = k*NC*n_neck;
n_long = length(k); % equal length of k
n_cont = 30;
n_run_ref = 1; % for runing random DOEs
mode = 6; % for call benchmark
h_len = 2;
lambda = NB;
%% Complement computation
Distance = zeros(n_neck, n_neck);
for i = 1:n_neck
    for j = i+1:n_neck
        Distance(i,j) = sqrt(string_kernel(Neck(i,:),Neck(i,:),Neck,lambda)+string_kernel(Neck(j,:),Neck(j,:),Neck,lambda)-2*string_kernel(Neck(i,:),Neck(j,:),Neck,lambda));
    end
end
Distance = Distance + (triu(Distance,1))';

% Compute the median d_neck
D = [];
for i = 1:n_neck
    tmp = [];
    for j = i+1:n_neck
        tmp = [tmp d_neck(Neck(i,:),Neck(j,:))];
    end
    D = [D tmp];
end
M_Y = median(D(:));
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
%% Target distribution (need X for each Neck)
n_tar = 10000;
F_ref = [];
for i = 1:n_neck
    X = zeros(n_tar, NC);
    L = rand(n_tar, NC);
    extent = upper - lower;
    for j=1:NC
        X(:,j) = extent(j)*L(:,j) + repmat(lower(j), n_tar, 1);
    end
    if mode == 1 || mode == 2 || mode == 3
        F_ref_tmp = zeros(NC,n_tar);
        for k = 1:n_tar
            F_ref_tmp(:,k) = kernel_benchmark_group(X(k,:),n_neck,i,mode);
        end
        F_ref = [ F_ref mean(F_ref_tmp,2)'];
    else
        F_ref_tmp = zeros(1,n_tar);
        for k = 1:n_tar
            F_ref_tmp(k) = kernel_benchmark_group(X(k,:),n_neck,i,mode);
        end
        F_ref = [F_ref mean(F_ref_tmp)];
    end
end
F_ref_mean = mean(F_ref);
eval(sprintf('save F_ref_mode_%d.mat F_ref;',mode))


Error_Rand_long = zeros(n_run_ref,n_long);
Error_LHS_long = zeros(n_run_ref,n_long);
Error_Sobol_long = zeros(n_run_ref,n_long);
Error_A_LHS_long = zeros(n_run_ref,n_long);
Error_MDS_long = zeros(n_run_ref,n_long);
Error_A_Greedy_long = zeros(n_run_ref,n_long);

F_Rand_save = [];
F_LHS_save = [];
F_A_LHS_save = [];
F_Sobol_save = [];
F_MDS_save = [];
F_A_Greedy_save = [];


%% Run R-Random and R-LHS 50 times for benchmarks
for i = 1:n_run_ref
    F_Rand_tmp = [];
    F_LHS_tmp = [];
    F_A_LHS_tmp = [];
    F_Sobol_tmp = [];
    F_MDS_tmp = [];
    F_A_Greedy_tmp = [];
    for h = 1:h_len:n_long
        X_Rand = zeros(nDOE(h), NC); % matrix load the values
        L = rand(nDOE(h), NC);
        extent = upper - lower;
        for j = 1:NC
            X_Rand(:,j) = extent(j)*L(:,j) + repmat(lower(j), nDOE(h), 1);
        end
        Y_Rand = round(rand(nDOE(h), NB));
        Pick_Rand = [Y_Rand X_Rand];
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
        % Sobol picking
        Pick_Sobol = LHS_cat(NC,NB,Neck,lower,upper,nDOE(h),0);
        % LHS picking
        Pick_A_LHS = LHS_cat(NC,NB,Neck,lower,upper,nDOE(h),1);
        % Greedy-MDS
        Pick_MDS = greedy_pick_cluster(NC,NB,Neck,lower,upper,nDOE(h),M_Y,Distance,K,lambda,n_cont,scale_x, 0);
        %greedy_pick_new(NC,NB,Neck,lower,upper,nDOE(h),M_Y,Distance, 0);
        % Adapted-Greedy
        Pick_A_Greedy = greedy_pick_cluster(NC,NB,Neck,lower,upper,nDOE(h),M_Y,Distance,K,lambda,n_cont,scale_x, 1);
        %greedy_pick_new(NC,NB,Neck,lower,upper,nDOE(h),M_Y,Distance,1);
        % Compute the subgroup and average of x in each method
        F_Rand = subgroup_and_compute_average_f(NC,NB,Pick_Rand,Neck,mode);
        F_LHS = subgroup_and_compute_average_f(NC,NB,Pick_LHS,Neck,mode);
        %%
        F_Sobol = subgroup_and_compute_average_f(NC,NB,Pick_Sobol,Neck,mode);
        F_A_LHS = subgroup_and_compute_average_f(NC,NB,Pick_A_LHS,Neck,mode);
        F_MDS = subgroup_and_compute_average_f(NC,NB,Pick_MDS,Neck,mode);
        F_A_Greedy = subgroup_and_compute_average_f(NC,NB,Pick_A_Greedy,Neck,mode);
        F_Rand_tmp = [F_Rand_tmp F_Rand];
        F_LHS_tmp = [F_LHS_tmp F_LHS];
        F_A_LHS_tmp = [F_A_LHS_tmp F_A_LHS];
        F_Sobol_tmp = [F_Sobol_tmp F_Sobol];
        F_MDS_tmp = [F_MDS_tmp F_MDS];
        F_A_Greedy_tmp = [F_A_Greedy_tmp F_A_Greedy];
       
        % Compute the errors: compute at each levels, then mean (or median errors of these)
        Err_Rand = 0;
        Err_LHS = 0;
        %%
        Err_Sobol = 0;
        Err_A_LHS = 0;
        Err_MDS = 0;
        Err_A_Greedy = 0;
        if mode == 1 || mode == 2 || mode == 3
            for k = 1:n_neck
                Err_Rand = Err_Rand + 1/NC*(norm(F_ref(:,k) - F_Rand(:,k)))^2;
                Err_A_LHS = Err_A_LHS + 1/NC*(norm(F_ref(:,k) - F_A_LHS(:,k)))^2;
                Err_LHS = Err_LHS + 1/NC*(norm(F_ref(:,k) - F_LHS(:,k)))^2;
                Err_Sobol = Err_Sobol + 1/NC*(norm(F_ref(:,k) - F_Sobol(:,k)))^2;
                Err_MDS = Err_MDS + 1/NC*(norm(F_ref(:,k) - F_MDS(:,k)))^2;
                Err_A_Greedy = Err_A_Greedy + 1/NC*(norm(F_ref(:,k) - F_A_Greedy(:,k)))^2;
            end
        else
            Err_Rand = (norm(F_ref - F_Rand))^2;
            Err_A_LHS =  (norm(F_ref - F_A_LHS))^2;
            Err_LHS = (norm(F_ref - F_LHS))^2;
            Err_Sobol =  (norm(F_ref - F_Sobol))^2;
            Err_MDS = (norm(F_ref - F_MDS))^2;
            Err_A_Greedy =  (norm(F_ref - F_A_Greedy))^2;
            
        end
        Error_Rand_long(i,h) = sqrt((1/n_neck)*Err_Rand);
        Error_LHS_long(i,h) =  sqrt((1/n_neck)*Err_LHS);
        %%
        Error_Sobol_long(i,h) =  sqrt((1/n_neck)*Err_Sobol);
        Error_A_LHS_long(i,h) =  sqrt((1/n_neck)*Err_A_LHS);
        Error_MDS_long(i,h) =  sqrt((1/n_neck)*Err_MDS);
        Error_A_Greedy_long(i,h) =  sqrt((1/n_neck)*Err_A_Greedy);
        
    end
    F_Rand_save = [F_Rand_save;F_Rand_tmp];
    F_LHS_save = [F_LHS_save;F_LHS_tmp];
    F_A_LHS_save = [F_A_LHS_save;F_A_LHS_tmp];
    F_Sobol_save = [F_Sobol_save;F_Sobol_tmp];
    F_MDS_save = [F_MDS_save;F_MDS_tmp];
    F_A_Greedy_save = [F_A_Greedy_save;F_A_Greedy_tmp];
    
end

Rand_mean = mean(Error_Rand_long);
LHS_mean = mean(Error_LHS_long);
Sobol_mean = mean(Error_Sobol_long);
A_LHS_mean = mean(Error_A_LHS_long);
MDS_mean = mean(Error_MDS_long);
A_Greedy_mean = mean(Error_A_Greedy_long);

Rand_median = median(Error_Rand_long);
LHS_median = median(Error_LHS_long);
Sobol_median = median(Error_Sobol_long);
A_LHS_median = median(Error_A_LHS_long);
MDS_median = median(Error_MDS_long);
A_Greedy_median = median(Error_A_Greedy_long);

Rand_std = std(Error_Rand_long);
LHS_std = std(Error_LHS_long);
Sobol_std = std(Error_Sobol_long);
A_LHS_std = std(Error_A_LHS_long);
MDS_std = std(Error_MDS_long);
A_Greedy_std = std(Error_A_Greedy_long);



%disp([ 1/n_neck*Err_Rand 1/n_neck*Err_LHS 1/n_neck*Err_Sobol 1/n_neck*Err_A_LHS 1/n_neck*Err_MDS  1/n_neck*Err_A_Greedy])
eval(sprintf(' save Rand_NB%d_NC%d_mode%d_soft.mat Error_Rand_long;',NB,NC,mode));
eval(sprintf(' save LHS_NB%d_NC%d_mode%d_soft.mat Error_LHS_long;',NB,NC,mode))
eval(sprintf(' save Sobol_NB%d_NC%d_mode%d_soft.mat Error_Sobol_long;',NB,NC,mode))
eval(sprintf(' save A_LHS_NB%d_NC%d_mode%d_soft.mat Error_A_LHS_long;',NB,NC,mode))
eval(sprintf(' save MDS_NB%d_NC%d_mode%d_soft.mat Error_MDS_long;',NB,NC,mode))
eval(sprintf(' save A_Greedy_NB%d_NC%d_mode%d_soft.mat Error_A_Greedy_long;',NB,NC,mode))

eval(sprintf(' save F_Rand_NB%d_NC%d_mode%d_soft.mat F_Rand_save;',NB,NC,mode));
eval(sprintf(' save F_LHS_NB%d_NC%d_mode%d_soft.mat F_LHS_save;',NB,NC,mode));
eval(sprintf(' save F_A_LHS_NB%d_NC%d_mode%d_soft.mat F_A_LHS_save;',NB,NC,mode));
eval(sprintf(' save F_Sobol_NB%d_NC%d_mode%d_soft.mat F_Sobol_save;',NB,NC,mode));
eval(sprintf(' save F_MDS_NB%d_NC%d_mode%d_soft.mat F_MDS_save;',NB,NC,mode));
eval(sprintf(' save F_A_Greedy_NB%d_NC%d_mode%d_soft.mat F_A_Greedy_save;',NB,NC,mode));

figure(1) 
errorbar(nDOE(1:h_len:n_long),Rand_mean(1:h_len:n_long),Rand_std(1:h_len:n_long),'m')
hold on
errorbar(nDOE(1:h_len:n_long),LHS_mean(1:h_len:n_long),LHS_std(1:h_len:n_long),'c')
hold on
errorbar(nDOE(1:h_len:n_long),Sobol_mean(1:h_len:n_long),Sobol_std(1:h_len:n_long),'k')
hold on
errorbar(nDOE(1:h_len:n_long),A_LHS_mean(1:h_len:n_long),A_LHS_std(1:h_len:n_long),'g')
hold on
errorbar(nDOE(1:h_len:n_long),MDS_mean(1:h_len:n_long),MDS_std(1:h_len:n_long),'b')
hold on
errorbar(nDOE(1:h_len:n_long),A_Greedy_mean(1:h_len:n_long),A_Greedy_std(1:h_len:n_long),'r')

legend('R-Random','R-LHS','A-Sobol','A-LHS','Greedy-MDS','Adapted-Greedy')



figure(2) 
errorbar(nDOE(1:h_len:n_long),Rand_median(1:h_len:n_long),Rand_std(1:h_len:n_long),'m')
hold on
errorbar(nDOE(1:h_len:n_long),LHS_median(1:h_len:n_long),LHS_std(1:h_len:n_long),'c')
hold on
errorbar(nDOE(1:h_len:n_long),Sobol_median(1:h_len:n_long),Sobol_std(1:h_len:n_long),'k')
hold on
errorbar(nDOE(1:h_len:n_long),A_LHS_median(1:h_len:n_long),A_LHS_std(1:h_len:n_long),'g')
hold on
errorbar(nDOE(1:h_len:n_long),MDS_median(1:h_len:n_long),MDS_std(1:h_len:n_long),'b')
hold on
errorbar(nDOE(1:h_len:n_long),A_Greedy_median(1:h_len:n_long),A_Greedy_std(1:h_len:n_long),'r')
legend('R-Random','R-LHS','A-Sobol','A-LHS','Greedy-MDS','Adapted-Greedy')



figure(3) % plot in logscale
errorbar(nDOE(1:h_len:n_long),log(Rand_mean(1:h_len:n_long)),(Rand_std(1:h_len:n_long)),'m')
hold on
errorbar(nDOE(1:h_len:n_long),log(LHS_mean(1:h_len:n_long)),(LHS_std(1:h_len:n_long)),'c')
hold on
errorbar(nDOE(1:h_len:n_long),log(Sobol_mean(1:h_len:n_long)),Sobol_std(1:h_len:n_long),'k')
hold on
errorbar(nDOE(1:h_len:n_long),log(A_LHS_mean(1:h_len:n_long)),A_LHS_std(1:h_len:n_long),'g')
hold on
errorbar(nDOE(1:h_len:n_long),log(MDS_mean(1:h_len:n_long)),MDS_std(1:h_len:n_long),'b')
hold on
errorbar(nDOE(1:h_len:n_long),log(A_Greedy_mean(1:h_len:n_long)),A_Greedy_std(1:h_len:n_long),'r')
legend('R-Random','R-LHS','A-Sobol','A-LHS','Greedy-MDS','Adapted-Greedy')

figure(4)
plot(nDOE(1:h_len:n_long),log(Rand_median(1:h_len:n_long)),'m')
hold on
plot(nDOE(1:h_len:n_long),log(LHS_median(1:h_len:n_long)),'c')
hold on
plot(nDOE(1:h_len:n_long),log(Sobol_median(1:h_len:n_long)),'k')
hold on
plot(nDOE(1:h_len:n_long),log(A_LHS_median(1:h_len:n_long)),'g')
hold on
plot(nDOE(1:h_len:n_long),log(MDS_median(1:h_len:n_long)),'b')
hold on
plot(nDOE(1:h_len:n_long),log(A_Greedy_median(1:h_len:n_long)),'r')
legend('R-Random','R-LHS','A-Sobol','A-LHS','Greedy-MDS','Adapted-Greedy')

%% run 1 times for adapted DOEs
% % Error_Sobol_long = zeros(n_long,1);
% % Error_A_LHS_long = zeros(n_long,1);
% % Error_MDS_long = zeros(n_long,1);
% % Error_A_Greedy_long = zeros(n_long,1);
% % 
% % for h = 1:n_long
% %     % Sobol picking
% %     Pick_Sobol = LHS_cat(NC,NB,Neck,lower,upper,nDOE(h),0);
% %     % LHS picking
% %     Pick_A_LHS = LHS_cat(NC,NB,Neck,lower,upper,nDOE(h),1);
% %     % Greedy-MDS
% %     Pick_MDS = greedy_pick(NC,NB,Neck,lower,upper,nDOE(h),M_Y,Distance, 0);
% %     % Adapted-Greedy
% %     Pick_A_Greedy = greedy_pick(NC,NB,Neck,lower,upper,nDOE(h),M_Y,Distance,1);
% %     
% %     % Compute the subgroup and average of x in each method
% %     X_mean_Sobol = subgroup_and_compute_average(NC,NB,Pick_Sobol,Neck);
% %     X_mean_A_LHS = subgroup_and_compute_average(NC,NB,Pick_A_LHS,Neck);
% %     X_mean_MDS = subgroup_and_compute_average(NC,NB,Pick_MDS,Neck);
% %     X_mean_A_Greedy = subgroup_and_compute_average(NC,NB,Pick_A_Greedy,Neck);
% %     
% %     % Approximation
% %     if mode == 1 || mode == 2 || mode == 3
% %         F_Sobol = zeros(n_neck*NC,1);
% %         F_A_LHS = zeros(n_neck*NC,1);
% %         F_MDS = zeros(n_neck*NC,1);
% %         F_A_Greedy = zeros(n_neck*NC,1);
% %         for i = 1:n_neck
% %             F_Sobol((i-1)*NC+1:i*NC) = kernel_benchmark(X_mean_Sobol(i,:),n_neck,i,mode);
% %             F_A_LHS((i-1)*NC+1:i*NC) = kernel_benchmark(X_mean_A_LHS(i,:),n_neck,i,mode);
% %             F_MDS((i-1)*NC+1:i*NC) = kernel_benchmark(X_mean_MDS(i,:),n_neck,i,mode);
% %             F_A_Greedy((i-1)*NC+1:i*NC) = kernel_benchmark(X_mean_A_Greedy(i,:),n_neck,i,mode);
% %         end
% %     else
% %         F_Sobol = zeros(n_neck,1);
% %         F_A_LHS = zeros(n_neck,1);
% %         F_MDS = zeros(n_neck,1);
% %         F_A_Greedy = zeros(n_neck,1);
% %         for i = 1:n_neck
% %             F_Sobol(i) = kernel_benchmark(X_mean_Sobol(i,:),n_neck,i,mode);
% %             F_A_LHS(i) = kernel_benchmark(X_mean_A_LHS(i,:),n_neck,i,mode);
% %             F_MDS(i) = kernel_benchmark(X_mean_MDS(i,:),n_neck,i,mode);
% %             F_A_Greedy(i) = kernel_benchmark(X_mean_A_Greedy(i,:),n_neck,i,mode);
% %         end
% %     end
% %     % Compute the errors: compute at each levels, then mean (or median errors of these)
% %     Err_Sobol = 0;
% %     Err_A_LHS = 0;
% %     Err_MDS = 0;
% %     Err_A_Greedy = 0;
% %     if mode == 1 || mode == 2 || mode == 3
% %         for i = 1:n_neck
% %             Err_Sobol = Err_Sobol + sqrt((1/NC*sum((F_ref((i-1)*NC+1:i*NC)-F_Sobol((i-1)*NC+1:i*NC)).^2)));
% %             Err_A_LHS = Err_A_LHS + sqrt((1/NC*sum((F_ref((i-1)*NC+1:i*NC)-F_A_LHS((i-1)*NC+1:i*NC)).^2)));
% %             Err_MDS = Err_MDS + sqrt((1/NC*sum((F_ref((i-1)*NC+1:i*NC)-F_MDS((i-1)*NC+1:i*NC)).^2)));
% %             Err_A_Greedy = Err_A_Greedy + sqrt((1/NC*sum((F_ref((i-1)*NC+1:i*NC)-F_A_Greedy((i-1)*NC+1:i*NC)).^2)));
% %         end
% %     else
% %         for i = 1:n_neck
% %             Err_Sobol = Err_Sobol + sqrt((1/NC*sum((F_ref(i)-F_Sobol(i)).^2)));
% %             Err_A_LHS = Err_A_LHS + sqrt((1/NC*sum((F_ref(i)-F_A_LHS(i)).^2)));
% %             Err_MDS = Err_MDS + sqrt((1/NC*sum((F_ref(i)-F_MDS(i)).^2)));
% %             Err_A_Greedy = Err_A_Greedy + sqrt((1/NC*sum((F_ref(i)-F_A_Greedy(i)).^2)));
% %         end
% %     end
% %         
% %     Error_Sobol_long(h) = (1/n_neck)*Err_Sobol;
% %     Error_A_LHS_long(h) = (1/n_neck)*Err_A_LHS;
% %     Error_MDS_long(h) = (1/n_neck)*Err_MDS;
% %     Error_A_Greedy_long(h) = (1/n_neck)*Err_A_Greedy;
% % end




