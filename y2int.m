function [XX] = y2int(x,y,argv)
% Rebuild original set of variables to be able to evaluate the function
% Converts binary variables into original integer variables with nlevels levels   
    XX = nan(argv.NC+argv.nZ,1);
    indexX = argv.indexX;
    indexZ = argv.indexZ;
    XX(indexX) = x;
    Nbin = argv.Nbin;   
    for i=1:argv.nZ
       Z = y(Nbin*(i-1)+1);
       for j=2:Nbin
          Z = Z + y(Nbin*(i-1)+j)*2^(j-1);
       end
       XX(indexZ(i)) = Z+1;
    end 
end