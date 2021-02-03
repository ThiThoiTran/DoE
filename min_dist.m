% SFD criteria: Using point-to-point and L_inf discrepancy
%% Point-to-point: using maximin
function m = min_dist(NC,NB,DOE)
[n_points,~] = size(DOE);
min_point = zeros(n_points,1);
for i = 1:n_points
    tmp = [];
    for  j = 1:n_points
        tmp = [tmp;d_neck(DOE(i,1:NB),DOE(j,1:NB))+ norm(DOE(i,NB+1:NB+NC)-DOE(j,NB+1:NB+NC))];
    end
    min_point(i) = min(tmp(find(tmp > eps)));
end
m = min(min_point);
end
        
    

