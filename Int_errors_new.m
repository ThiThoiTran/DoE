%% Run different types of DOEs in a set of benchmark functions 
% With new measurament: expectation
close all
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
MC_ref = dlmread('MCIE_Vegas_2.txt');
M_dim = [2 2;4 3;6 2;10 4;20  6;1 3;2 3; 5 3;5 3;7 3;14 3;2 8;7 8;1 8;1 10];
M_bound = repmat([0 2],n_prob,1);
M_level = [3 4 3 6 14 4 4 4 4 4 4 36 36 36 108] ;
%% Create and compute the error for each function
mode = 2;
name = eval(sprintf('name_%d',mode));
NC = eval(sprintf('M_dim(%d,1)',mode));
NB = eval(sprintf('M_dim(%d,2)',mode));
dim = NB + NC;
Neck = countSol(NB);
n_neck = size(Neck,1);
n_repeat = 1;
k = 1:50:500;
nDOE =  k*(NC+size(Neck,1));%[10 100:100:1000]; %
n_test = length(nDOE);
lower_x = eval(sprintf('M_bound(%d,1)',mode));
upper_x = eval(sprintf('M_bound(%d,2)',mode));
lower = [lower_x*ones(1,NC) zeros(1,NB)];
upper = [upper_x*ones(1,NC) ones(1,NB)];

%% Build the empirical distribution
if NB <= 10
    n_cont = 100;
else
    n_cont = 30;
end
X = zeros(n_cont*n_neck, NC);
for i = 1:n_neck
    L = lhsdesign(n_cont, NC);
    extent = upper - lower;
    for j=1:NC
        X((i-1)*n_cont+1:i*n_cont,j) = extent(j)*L(:,j) + repmat(lower(j), n_cont, 1);
    end
end
Mat = [kron(Neck, ones(n_cont,1)) X];
f_ref = OF(NC,NB,Mat,Neck,mode);

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
        tmp = OF(NC,NB,[Y((k-1)*nDOE(i)+1:k*nDOE(i),:) X((k-1)*nDOE(i)+1:k*nDOE(i),:)],Neck,mode);
        err_tmp = abs(f_ref-tmp)/abs(f_ref);
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
        tmp = OF(NC,NB,[Y((k-1)*nDOE(i)+1:k*nDOE(i),:) X((k-1)*nDOE(i)+1:k*nDOE(i),:)],Neck,mode);
        err_tmp = abs(f_ref-tmp)/abs(f_ref);
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
        tmp = OF(NC,NB,[Y((k-1)*nDOE(i)+1:k*nDOE(i),:) X((k-1)*nDOE(i)+1:k*nDOE(i),:)],Neck,mode);
        err_tmp = abs(f_ref-tmp)/abs(f_ref);
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
        tmp = OF(NC,NB,[Y((k-1)*nDOE(i)+1:k*nDOE(i),:) X((k-1)*nDOE(i)+1:k*nDOE(i),:)],Neck,mode);
        err_tmp = abs(f_ref-tmp)/abs(f_ref);
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
        tmp = OF(NC,NB,[Y((k-1)*nDOE(i)+1:k*nDOE(i),:) X((k-1)*nDOE(i)+1:k*nDOE(i),:)],Neck,mode);
        err_tmp = abs(f_ref-tmp)/abs(f_ref);
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
%%



%% 
% Save Error results
% save Error10repets.mat Rand_error_simulation Rand_Er_simulation LHS_error_simulation LHS_Er_simulation Sobol_error_simulation Sobol_Er_simulation MDS_error_simulation MDS_Er_simulation Adapt_error_simulation Adapt_Er_simulation
fprintf('\n R-Random  \t  R-LHS \t  A-Sobol \t A-LHS \t  MDS \n');
Error = [eval(sprintf('Rand_error_%s(:)',eval(sprintf('name_%d',mode)))) eval(sprintf('LHS_error_%s(:)',eval(sprintf('name_%d',mode))))...
    eval(sprintf('Sobol_error_%s(:)',eval(sprintf('name_%d',mode)))) eval(sprintf('A_LHS_error_%s(:)',eval(sprintf('name_%d',mode))))...
    eval(sprintf('MDS_error_%s(:)',eval(sprintf('name_%d',mode))))];%  eval(sprintf('Adapt_error_%s(:)',eval(sprintf('name_%d',mode))))];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(1)
% subplot(1,3,1)
% bar(ind_rand(:));
% title('Random DOEs')
% subplot(1,3,2)
% bar(ind_lhs(:));
% title('LHS DOEs')
% subplot(1,3,3)
% bar(ind_A_LHS(:));
% title('Adapted-LHS DOEs')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(2)
% subplot(1,3,1)
% bar(ind_sobol)
% title('Sobol sequence DOEs')
% subplot(1,3,2)
% bar(ind_mds)
% title('Greedy-MDS DOEs')
% subplot(1,3,3)
% bar(ind_A_LHS)
% title('A-LHS DOEs')

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
% scatter(6*ones(n_repeat,1),eval(sprintf('Adapt_error_%s(:)',eval(sprintf('name_%d',mode)))),'g','filled','jitter','on','jitterAmount',0.1)
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
plot(nDOE,eval(sprintf('Rand_Er_%s(:)',eval(sprintf('name_%d',mode)))),'c')
hold on
plot(nDOE,eval(sprintf('LHS_Er_%s(:)',eval(sprintf('name_%d',mode)))),'g')
hold on
plot(nDOE,eval(sprintf('Sobol_Er_%s(:)',eval(sprintf('name_%d',mode)))),'k')
hold on
plot(nDOE,eval(sprintf('MDS_Er_%s(:)',eval(sprintf('name_%d',mode)))),'b')
% hold on
% plot(nDOE,eval(sprintf('Adapt_Er_%s(:)',eval(sprintf('name_%d',mode)))),'r')
hold on
plot(nDOE,eval(sprintf('A_LHS_Er_%s(:)',eval(sprintf('name_%d',mode)))),'r')
legend('R-Random','R-LHS','A-Sobol','Greedy-MDS','A-LHS')
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