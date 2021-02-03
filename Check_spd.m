% Check if a matrix is spd:
NB_max = 14;
NB_min = 2;
n = NB_max - NB_min + 1;
lambda = 0.1;
for i = 6 % 6:NB_max
    NB = i;
    Neck = countSol(i);
    n_neck = size(Neck,1);
    D = zeros(n_neck,n_neck);
    D1 = zeros(n_neck,n_neck);
    for j = 1:n_neck
        for k = j+1:n_neck
            D(j,k) = exp(-lambda*(d_neck(Neck(j,:),Neck(k,:)))-log(n_neck));
            D1(j,k) = exp(-sum(abs(Neck(j,:)-Neck(k,:))));
        end
    end
    D = D + D'+ eye(n_neck);
    D1 = D1 + D1'+ eye(n_neck);
    disp([min(eig(D)) min(eig(D1))]);
    if min(eig(D))>=0
        fprintf('Matrix size %d is SDP\n',i);
    else
        fprintf('Matrix size %d is not SDP\n',i);
    end
end




    
            
    





