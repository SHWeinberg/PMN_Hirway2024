%% CPM Trial 4.2- Wound Healing Assay
% Run confluent monolayer with scratch testing- cells in the middle
% disappear

%% Load Simulation
% get ctag matrix and intra/extracellular concentrations
load('test_100p_1201mcs_5tgfb_0.3tmax_0.06spread_targetscale2_Jhalf9574.3821_dividing1.25_extraA0.96_extraBT0.6_ZEBscale0.8_R200scale1.2_ZEB1_SNAIL1_TGFB.mat') %_inelasticity_1000
finalctag=ctag(:,end);
initialctag=ctag(:,1);
dimval=100;
finalgridconc=statevars{end}(1:3*dimval*dimval);
% finalctag=reshape(finalctag,100,100);
%% Scratch Test
finaltag2=reshape(finalctag,dimval,dimval);
center=dimval/2;
cut=15;
figure
subplot(1,2,1)
imagesc(finaltag2);

title('Before');
finaltag3=finaltag2;
finaltag3(center-cut:center+cut-1,:)=0; % makes the cut, turning cells to 0

NRc=max(max(finaltag2));
csize=zeros(NRc,3);
for cell=1:NRc
    for c = 1:dimval*dimval
        if (finaltag2(c) == cell)
            csize(cell,1) = csize(cell,1) + 1; % measures cell size before the cut
        end
    end
end

for cell=1:NRc
    for c = 1:dimval*dimval
        if (finaltag3(c) == cell)
            csize(cell,2) = csize(cell,2) + 1; % measures cell size after the cut
        end
    end
end
csize(:,3)=csize(:,1)-csize(:,2); % finds difference in before and after

% Reorganize ctag & grid
eliminated=(csize(:,1)-csize(:,2)>0); %Removes cells that are even partially cut
elimpos=find(eliminated==1);
finaltag4=reshape(finaltag3,dimval*dimval,1);
% Redo the grid to remove partial cells
for e=1:size(elimpos,1)
    for i =1: size(finaltag4,1)
        if (finaltag4(i)==elimpos(e))
            finaltag4(i)=0;
        end
        
    end
end

disp('Partials cut');
% reorder ctag to remove gaps in cell # order
for e=size(elimpos,1):-1:1
    for i =1: size(finaltag4,1)
        if (finaltag4(i)>=elimpos(e))
            finaltag4(i)=finaltag4(i)-1;
        end
        
    end
end
NRc=max(max(finaltag4));
csize2=zeros(NRc,3);
for cell=1:NRc
    for c = 1:dimval*dimval
        if (finaltag4(c) == cell)
            csize2(cell,1) = csize2(cell,1) + 1;
        end
    end
end
subplot(1,2,2)
imagesc(reshape(finaltag4,dimval,dimval));
title('After');

% Reorder state variables
stateV=statevars{end}(3*dimval*dimval+1:end);
stateV2=reshape(stateV,size(stateV,1)/8,8);
stateV2(elimpos,:)=[];

%% Set up Initial conditions
dimval=100;
mcsval=150;
tgfbval=0;
spreadval=0.015*4;
tmaxval=0;
Jmin=4819*15;
dividing=1.25;%/4;
al=0.05;
n=2;
Jhalf0=Jmin/((al^-1-1)^(1/n));
Tmax=2*.15;
tmax0=1*.15;
Jmin= Jhalf0*(al^-1 -1) ^(1/n) / (2*Tmax/tmax0 - 1)^(1/n);
Jhalf=Jmin/((al^-1-1)^(1/n));

pZEB=0.8;
pR200=1.2;
tmaxval=Tmax;
extrascale1=0.96;
extrascale2=0.6; % ECM TGFB Binding Rate

intrazeb1=1;
intrasnail2=1;
kdae_scale=1;
comment=['Mcell&Egrid&scratch' num2str(cut)];


% Initial Seeding Pattern

% Use ctag from a simulation
initmatrix=finaltag4;

%Initial Seeding cell concentrations
icstate = [.01; .01; .38; ...
    .03; .01; .35; 3.2; 0];  % initial state, endoTGFB 0.16 omitted, NEED TO KEEP THIS
%cellvars=[];
% for c=1:max(max(initmatrix))
%     cellvars=[cellvars, icstate];
% end
% cellvars=cellvars';

% Use cell final statevars from simulation (Could use epithelial or
% mesenchymal here)
cellvars=stateV2;
cellvars=reshape(cellvars,max(max(initmatrix))*size(icstate,1),1);

% Initial seeding matrix concentrations
imatstate=[0.1571,0.0072,0.01,2.4569*10^-5];
% matrvars=nan(dimval*dimval,4);
% matrvars(:,1)=imatstate(1);
% matrvars(:,2)=imatstate(2);
% matrvars(:,3)=imatstate(3);
% matrvars(:,4)=imatstate(4);
% matrvars=reshape(matrvars,4*dimval*dimval,1);

griddim=dimval*dimval;
% Alternatively, get grid concentrations from simulation, (Could use
% epithelial or mesenchymal grid here)
matrvars=nan(dimval*dimval,3);
matrvars(:,1)=finalgridconc(1:1*dimval*dimval);
matrvars(:,2)=finalgridconc(1*dimval*dimval+1:2*dimval*dimval);
matrvars(:,3)=finalgridconc(2*dimval*dimval+1:3*dimval*dimval);
matrvars=reshape(matrvars,3*dimval*dimval,1);

%% Change grid conc for scratch
for i=1:10000
    if(initmatrix(i)==0)
        matrvars(i)=epigridconc(i);
        matrvars(10000+i)=epigridconc(10000+i);
        matrvars(20000+i)=epigridconc(0000+i);
    end
end
%%
clc

run_single_cpm_fem_ode(dimval,mcsval,tgfbval,spreadval,tmaxval,Jhalf,dividing,extrascale1,extrascale2,intrazeb1, intrasnail2,pZEB,pR200,kdae_scale,comment,initmatrix,cellvars,matrvars);
pause(1)
close all

