% Create initital sample set based on LHSmaximin in the necklace space
function S = LHS_cat(NC,Nblades,levels,lower,upper,hnum,type)
% NC: number of continuous variables
% Nblades: Number of binary variables originally
% n_sample: number of set to create
% name: name to save in .mat file

NB = Nblades;
argv.lower = lower;
argv.upper = upper;
X = zeros(hnum, NC); % matrix load the values
n_multi = floor(hnum/size(levels,1));
n_add = hnum - floor(hnum/size(levels,1))*size(levels,1);  

%% Initial points from
weights = 1/size(levels,1);
L = lhsdesign(hnum, NC);
extent = argv.upper - argv.lower;
for j=1:NC
    X(1:hnum,j) = extent(j)*L(:,j) + repmat(argv.lower(j), hnum, 1);
end
Y = [];
if ~n_multi
    for i = 1:n_multi
        if type == 0 % Using sobol sequence to create points
            p = sobolset(1,'Skip',1e3,'Leap',1e2);
            p = scramble(p,'MatousekAffineOwen');
            var_cat = net(p,size(levels,1));
        elseif type == 1 % Using LHS to create points
            var_cat = lhsdesign(size(levels,1), 1);
        end
        vec_cat = con_to_cat(var_cat,levels,weights);
        temp = cat_to_bi(vec_cat,levels,NB);
        Y = [Y;temp];
    end

    if n_add ~= 0
        if type == 0 % Using sobol sequence to create points
            p = sobolset(1,'Skip',1e3,'Leap',1e2);
            p = scramble(p,'MatousekAffineOwen');
            var_cat = net(p,n_add);
        elseif type == 1 % Using LHS to create points
            var_cat = lhsdesign(n_add, 1);
        end
        vec_cat = con_to_cat(var_cat,levels,weights);
        temp_1 = cat_to_bi(vec_cat,levels,NB);
        Y = [Y;temp_1];
    end
else
    if type == 0 % Using sobol sequence to create points
        p = sobolset(1,'Skip',1e3,'Leap',1e2);
        p = scramble(p,'MatousekAffineOwen');
        var_cat = net(p,hnum);
    elseif type == 1 % Using LHS to create points
        var_cat = lhsdesign(hnum, 1);
    end
    vec_cat = con_to_cat(var_cat,levels,weights);
    Y = cat_to_bi(vec_cat,levels,NB);
end
S = [ Y X];
end
%
% %% Function to convert contiuous to categorical variable
function [s] = con_to_cat(var_cat,levels,weights)
n_levels = size(levels,1);
if length(weights)==1
    weights = repmat(weights,n_levels,1);
end
n_points = length(var_cat);
interval = [0; cumsum(weights)];
s = zeros(n_points,1);
for i = 1:n_levels-1
    ind1 = find( var_cat <= interval(i+1));
    ind2 = find( var_cat > interval(i));
    ind = intersect(ind1,ind2);
    s(ind) = i;
end
ind = find(var_cat > interval(n_levels));
s(ind) = n_levels;
end
%% function to convert categorical variable to binary
function [y] = cat_to_bi(vec_cat,levels, NB)
n_cat = length(vec_cat);
for i = 1:n_cat
    y(i,:) = levels(vec_cat(i),:);
end
end

%% function to convert cate to necklace binary
function [y] = cat_to_neck(vec_cat,n)
n_cat = length(vec_cat);
y = zeros(n_cat,n);
for i = 1:n_cat
    [~,e]=log2(vec_cat(i));
    T =  rem(floor(vec_cat(i)*pow2(1-max(n,e):0)),2);
    y(i,:) = flip(T);
end
end




