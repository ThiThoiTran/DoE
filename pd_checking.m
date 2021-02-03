% Numeriacally check the positive definite of the necklace distance
% Idea of checking: Generate input: x,y, compute the matrix of distance,
% hence generate randomly vectors c, check if c^T K c >= 0
function pd_checking(NC,NB,n_check, n_c)
%% Generating input: (note that the input are vectors of binary)
% Choose the dimension
%NC = 1;
%NB = 6;
nDOE = NB+NC+1; 
% Choose number of total input checking
%n_check = 100;
Y = zeros(nDOE*n_check, NB);
for i = 1:n_check
    Y((i-1)*nDOE +1:i*nDOE,:) = randi([0 1], nDOE, NB);
end
%% Checking if the kernel with the necklace distance is positve of not
% Create randomly checking vectors c \in R^nDOE
%n_c = 100; % for each nDOE input, check 100 c
checking = zeros(n_c*n_check,1);
c_lower = -100;
c_upper = 100;
for i = 1: n_check
    K = zeros(nDOE, nDOE);
    Y_distract = Y((i-1)*nDOE+1:i*nDOE,:);
    for j = 1:nDOE
        for k = i+1:nDOE
            K(i,j) = exp(-1/2*sum(d_neck(Y_distract(i,:), Y_distract(j,:))));
          %   K(i,j) = exp(-1/2*sum(abs(Y_distract(i,:)-Y_distract(j,:)))); 
        end
    end 
    K = K + K';
    rng('shuffle')
    for h = 1:n_c
        c = (c_upper - c_lower)*rand(nDOE,1) + c_lower;
        t = c'*K*c;
        if t>=0
            checking((i-1)*n_c+h) = 1;
        end
    end
end
index = find(checking >= 0);
l = length(index);
if l==n_c*n_check
    fprintf('The kernel is numerically positive definite')
else
    fprintf('The kernel is not positive definite')
end



