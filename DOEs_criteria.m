
%% Checking maxiimn criteria: running 100 repetitions to see
close all
NB = 6;
NC = 2;
NB_tune = 2;
Neck = countSol(NB);
x_l = 0;
x_u = 5;
lower = [x_l*ones(1,NC) zeros(1,NB)];
upper = [x_u*ones(1,NC) ones(1,NB)];
[n_neck, n_dim] = size(Neck);
k = 5;
nDOE = k*(NC+n_neck);
n_cont = 30;
n_take = n_neck;
n_repeat = 30;
mindist_Rand = zeros(n_repeat,1);
mindist_R_LHS = zeros(n_repeat,1);
mindist_Sobol = zeros(n_repeat,1);
mindist_LHS = zeros(n_repeat,1);
mindist_MDS = zeros(n_repeat,1);
mindist_Adapted = zeros(n_repeat,1);
for niter = 1:n_repeat
    L = rand(nDOE, NC);
    extent = upper - lower;
    X_rand = zeros(nDOE,NC);
    for j=1:NC
        X_rand(:,j) = extent(j)*L(:,j) + repmat(lower(j), nDOE, 1);
    end
    Y_rand = round(rand(nDOE, NB));
    Rand_pick = [Y_rand X_rand];
    
    X_LHS = zeros(nDOE,NC);
    L1 = lhsdesign(nDOE, NC);
    extent = upper - lower;
    for j=1:NC
        X_LHS(:,j) = extent(j)*L1(:,j) + repmat(lower(j), nDOE, 1);
    end
    Y_LHS = round(lhsdesign(nDOE, NB));
    R_LHS_pick = [Y_LHS X_LHS];
    % Sobol picking
    Sobol_pick = LHS_cat(NC,NB,Neck,lower,upper,nDOE,'.',0);
    % LHS picking
    A_LHS_pick = LHS_cat(NC,NB,Neck,lower,upper,nDOE,'.',1);
    % Greedy-MDS picking
    pick_MDS = greedy_pick(NC,NB,Neck,lower,upper,nDOE,0);
    % Adapt mixed kernel directly to Greedy algorithm
    [pick_adapt] = greedy_pick(NC,NB,Neck,lower,upper,nDOE,1);
    
    mindist_Rand(niter) = min_dist(NC,NB,Rand_pick);
    mindist_R_LHS(niter) = min_dist(NC,NB,R_LHS_pick);
    mindist_Sobol(niter) = min_dist(NC,NB,Sobol_pick);
    mindist_LHS(niter) = min_dist(NC,NB,Sobol_pick);
    mindist_MDS(niter) = min_dist(NC,NB,pick_MDS);
    mindist_Adapted(niter) = min_dist(NC,NB,pick_adapt);
end
mindist = [mindist_Sobol mindist_LHS  mindist_MDS mindist_Rand mindist_R_LHS];
figure(1)
scatter(ones(n_repeat,1),mindist_Sobol,'r','filled','jitter','on','jitterAmount',0.1)
hold on
scatter(2*ones(n_repeat,1),mindist_LHS,'b','filled','jitter','on','jitterAmount',0.1)
hold on
scatter(3*ones(n_repeat,1),mindist_MDS,'k','filled','jitter','on','jitterAmount',0.1)
% hold on
% scatter(4*ones(n_repeat,1),mindist_Adapted,'y','filled','jitter','on','jitterAmount',0.1)
hold on
scatter(4*ones(n_repeat,1),mindist_Rand,'c','filled','jitter','on','jitterAmount',0.1)
hold on
scatter(5*ones(n_repeat,1),mindist_R_LHS,'g','filled','jitter','on','jitterAmount',0.1)

boxplot(mindist);
names = {'Adapted-Sobol';  'Adapted-LHS'; 'Greedy-MDS';'Rounding Random';'Rounding LHS'};
set(gca,'xtick',[1:5],'xticklabel',names,'Fontsize',8,'LineWidth',1)
ylabel('maximin')
title('Maximin criteria')
