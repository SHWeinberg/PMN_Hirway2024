function [J]=Jccpheno(tag1,tag2,phenotype,phenoparams,cellmod)
Jcctable=cellmod.JccTab;
cellAtype=cellmod.pop_index(tag1);
cellBtype=cellmod.pop_index(tag2);

cellA=phenotype(tag1);
cellB=phenotype(tag2);
n=4;
phalf=0.5;
E=phenoparams.jccE;
M=phenoparams.jccM;
cellAphen=0;
cellBphen=0;
% A_mostM = E + (M-E)*max(P1,P2).^n./(max(P1,P2).^n + phalf^n);
% A_mostE = E + (M-E)*min(P1,P2).^n./(min(P1,P2).^n + phalf^n);
%A_prod2 = E + (M-E)*(P1+P2-P1.*P2).^n./((P1+P2-P1.*P2).^n + phalf^n);
%A_prod = E + (M-E)*(P1.*P2).^n./((P1.*P2).^n + phalf^n);
% J=A_mostM;

%Cell 1
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


% E-E
% % if(0<=cellA && cellA<=0.48 && 0<=cellB && cellB<=0.48)
% %     J=1.75;
% % %disp('Case 1')
% % %E-P,
% % elseif((0<=cellA && cellA<=0.48 && 0.48<cellB && cellB<0.76) || (0.48<cellA && cellA<0.76 && 0<=cellB && cellB<=0.48))
% %     J=1.75;
% %  %   disp('Case 2')
% % %E-M,
% % elseif((0<=cellA && cellA<=0.48 && 0.76<=cellB) || (0.76<=cellA && 0<=cellB && cellB<=0.48))
% %     J=1.75;
% %   %  disp('Case 3')
% % %P-P,
% % elseif(0.48<cellA && cellA<0.76 && 0.48<cellB && cellB<0.76)
% %     J=1.425;
% %    % disp('Case 4')
% % %P-M,
% % elseif((0.48<cellA && cellA<0.76 && 0.76<=cellB) || (0.76<=cellA && 0.48<cellB && cellB<0.76))
% %     J=1.425;
% %     %disp('Case 5')
% % %M-M
% % else %if(2/3<=cellA && cellA<=1 && 2/3<=cellB && cellB<=1)
% %     J=1.1;
% %     %disp('Case 6')
% % end

J=Jcctable{cellAtype,cellBtype}(cellAphen,cellBphen);

end