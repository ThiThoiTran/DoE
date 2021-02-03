% function compute mean of function at a given DOEs
function [m,index1] = mean_DOE(NC,NB,DOEs,Ref,nDOE,mode,nlevel)
% Need to compute each level
index = zeros(nlevel,1);
F_sum = zeros(nlevel,1);
for i = 1:nDOE
    for  j = 1:nlevel 
            if myisrotation(DOEs(i,1:NB),Ref(j,:))
                index(j) = index(j) + 1;
                F_sum(j) = F_sum(j) + benchmark(mode,DOEs(i,NB+1:NB+NC),DOEs(i,1:NB),Ref);
                break;
            end
    end
end
n_visited_neck = numel(find(index>0));
index1 = index;
if all(index)
    m = F_sum./(index);
else
   index(find(index == 0)) = 1;
    m = F_sum./(index);
end
    
%m = sum(m)/nlevel;%n_visited_neck;
end
   