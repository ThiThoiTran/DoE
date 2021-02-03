close all
NB = 6;
NC = 1;
NB_tune = 2;
Neck = countSol(NB);
lower = [0*ones(1,NC) zeros(1,NB)];
upper = [1*ones(1,NC) ones(1,NB)];
[n_neck, n_dim] = size(Neck);
k = 1;
nDOE = k*(NC+n_neck);%NB +NC +1;
n_cont = 30;
n_take = n_neck;
X = zeros(n_cont*n_take, NC);
for i = 1:n_take
    L = lhsdesign(n_cont, NC);
    extent = upper - lower;
    for j=1:NC
        X((i-1)*n_cont+1:i*n_cont,j) = extent(j)*L(:,j) + repmat(lower(j), n_cont, 1);
    end
end

Mat1 = [kron(Neck, ones(n_cont,1)) X];
[pick_adapt] = Greedy(NC,NB,Mat1,nDOE,1);
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
figure(1)
scatterhist(pick_adapt(:,end),cat_adapt)
figure(2)
scatter(Mat1(:,end),cat(:),'b')
hold on
scatter(pick_adapt(:,end),cat_adapt,'r')

num_cat_adapt = zeros(n_neck,1);
for i = 1:n_neck
    num_cat_adapt(i) = numel(find(cat_adapt == i));
end

figure(3)
bar(num_cat_adapt)


