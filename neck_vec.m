% Create the vector of dependent location of points
function [v,n_cont] = neck_vec(n_neck,n_cont)
 tmp = n_neck - 2;
if rem(tmp/2,2) == 0
    v = [ repmat([1 2],1,tmp/4 ) 1];
else
    v = [ repmat([1 2],1,floor(tmp/4)) 2];
end
v = [v(end:-1:1) v ];

%v =  [n_cont*[n_neck/2:-1:2] 5 5 n_cont*[2:1:n_neck/2]]; %n_cont*[n_neck/2:-1:1 1:1:n_neck/2];
 v = n_cont*v;
% n_cont = sum(v);
n_cont = sum(v);
 
end