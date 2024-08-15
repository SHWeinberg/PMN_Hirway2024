function [J]=Jcmpheno(tag1,phenotype,phenoparams,cellmod)
Jcmtable=cellmod.JcmTab;
cellAtype=cellmod.pop_index(tag1);
cellA=phenotype(tag1);
n=2;
phalf=0.5;
E=phenoparams.jcmE;
M=phenoparams.jcmM;
%  J = E + (M-E)*(P1).^n./((P1).^n + phalf^n);
cellAphen=0;
%E 
if(0 <= cellA && cellA< 0.48)
    cellAphen=1;
   % J=2;
%disp('Case 1')
%P
elseif(0.48 <= cellA && cellA< 0.76)
    cellAphen=2;
    %J=1.5;
 %   disp('Case 2')
%M
else %if(2/3 <= cellA && cellA <=1)
    cellAphen=3;
   % J=1;
  %  disp('Case 3')

end
J=Jcmtable{cellAtype}(cellAphen);
end