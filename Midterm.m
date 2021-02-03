% File rum DFOb with dnew, with several initial set of points
close all
%%  Run 300 times SAFRAN full simulation 
% n_take = 100;
% for i = 1:n_take
%     fprintf('###################################### Sample set number %d',i)
%     eval(sprintf('Default_run(''620'',%d,''.'');',i));
%     load('store_620.mat','res');
%     eval(sprintf('save store_620_simul_full_sobol_ini_3_%d.mat res;',i));
% end

 n_probs = 25; % 25 benchmark functions
 n_take = 10; % each functions run for 10 repeatations
 name_1 = 'CB2';
 name_2 = 'CB3';
 name_3 = 'QL';
 name_4 = 'WF';
% name_5 = 'EVD52';
 name_5 = 'Branin';
 name_6 = 'pentagon';
 name_7 = 'rosensuzuki';
 name_8 = 'wong2';
 name_9 = 'wong3';
 name_10 = 'MAD1';
 name_11 = 'MAD2';
 %name_12 = 'HS78';
 name_12 = 'sporttournament';
 name_13 = 'HS2';
 name_14 = 'HS3';
 name_15 = 'HS29log';
name_16 = 'Hartman3';
name_17 = 'Hartman6';
name_18 = 'Sixhump';
name_19 = 'Goldstein';
%name_6 = 'Shekel5';
name_20 = 'Shekel7';
name_21 = 'Shekel10';
name_22 = 'ex811';
name_23 = 'ex814';
%name_11 = 'least';
name_24 = 'perm6';
name_25 = 'perm8';
%name_14 = 'miqp1';
%name_15 = 'sporttournament';
 
 for i = 19:n_probs
     %   eval(sprintf('Default_run(''620_%s'',1,''.'');',eval(sprintf('name_%d',i))));
     for j = 1:n_take
        eval(sprintf('Default_run(''620'',j,name_%d);',i));
        
        load('store_620.mat','res');
        eval(sprintf('save Sobol_%s_%d.mat res;',eval(sprintf('name_%d',i)),j));
%         SV = res.storeF;
%         [val, ind] = min(SV);
%         Imin =  ind;
%         Store = [ res.storeX(ind,:) res.storeY(ind,:) res.storeF(ind)] ;
%         M = [];
%         for j = 1:size(SV,1)
%             M = [M  min(SV(1:j))];
%         end
%         %V = [V V1];
%         n_simtot =  res.argv.simtot;
%         eval(sprintf('save Test_DOE_d_neck_%s.mat M Imin Store  n_simtot;',eval(sprintf('name_%d',i))));
     end
 end
   
    

 
 
 
