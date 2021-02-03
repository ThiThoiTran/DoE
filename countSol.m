function [ListofSol] = countSol(N)
% count the number of feasible disctinct solutions of 2 different shapes of
% blades among N blades

ListofSol = permn(0:1,N);
Nlist = size(ListofSol,1);

keep=ones(Nlist,1);
for i=1:Nlist-1
    ref = ListofSol(i,:);
    for j=i+find(keep((i+1):Nlist)==1)'
        k=1;
        while (keep(i) && k < N)
            jk = [ (1:N-k)+k 1:k ]; 
            ind = norm(ref - ListofSol(j,jk),1);
            keep(j) = (~(ind==0)) * keep(j);
            k=k+1;
        end
    end
end

temp=ListofSol;
ListofSol=[];
for i=1:Nlist
    if keep(i)
        ListofSol = [ListofSol ; temp(i,:)];
    end
end

%fprintf('%d Solutions \n',size(ListofSol,1));

end