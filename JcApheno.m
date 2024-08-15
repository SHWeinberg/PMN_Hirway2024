function [J]=JcApheno(tag1,tag2,phenotype,phenoparams,cellmod)
Jcatable=cellmod.JcaTab;
cellAtype=cellmod.pop_index(tag1);
cellBtype=cellmod.pop_index(tag2);

cellA=phenotype(tag1);
cellB=phenotype(tag2);

cellAphen=0;
cellBphen=0;

if(0 <= cellA && cellA< 0.48)
    cellAphen=1;
elseif(0.48 <= cellA && cellA< 0.76)
    cellAphen=2;
else %if(2/3 <= cellA && cellA <=1)
    cellAphen=3;
end

%Cell 2
if(0 <= cellB && cellB< 0.48)
    cellBphen=1;
elseif(0.48 <= cellB && cellB< 0.76)
    cellBphen=2;
else %if(2/3 <= cellA && cellA <=1)
    cellBphen=3;
end
%n=2;
%phalf=0.5;
%E=phenoparams.jccE;
%M=phenoparams.jccM;
% J = E + (M-E)*max(P1,P2).^n./(max(P1,P2).^n + phalf^n);
%disp('Jca')
%E-E, 
% % if(0<=cellA && cellA<=0.48 && 0<=cellB && cellB<=0.48)
% %     J=0.5;
% % %disp('Case 1')
% % %E-P,
% % elseif((0<=cellA && cellA<=0.48 && 0.48<cellB && cellB<0.76) || (0.48<cellA && cellA<0.76 && 0<=cellB && cellB<=0.48))
% %     J=0.3;
% %   %  disp('Case 2')
% % %E-M,
% % elseif((0<=cellA && cellA<=0.48 && 0.76<=cellB) || (0.76<=cellA && 0<=cellB && cellB<=0.48))
% %     J=0.1;
% %   %  disp('Case 3')
% % %P-P, 
% % elseif(0.48<cellA && cellA<2/3 && 0.48<cellB && cellB<0.76)
% %     J=0.175;
% %  %   disp('Case 4')
% % %P-M, 
% % elseif((0.48<cellA && cellA<0.76 && 0.76<=cellB) || (0.76<=cellA && 0.48<cellB && cellB<0.76))
% %     J=0.05;
% %   %  disp('Case 5')
% % %M-M
% % else %if(2/3<=cellA && cellA<=1 && 2/3<=cellB && cellB<=1)
% %     J=0;
% %  %   disp('Case 6')
% % end

%Removing Jca component SUH 1/26/20
%J=0; JcA restored SUH 1/28/20


J=Jcatable{cellAtype,cellBtype}(cellAphen,cellBphen);
end