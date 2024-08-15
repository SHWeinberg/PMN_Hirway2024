
Emat='test_100p_1200mcs_0.1tgfb_0.06spread1.25div_extraA0.96_extraBT0.6_ZEB1.125_SNAIL1_kDAE.scale1_Egrid.mat';
load(Emat)

%%
dimval=100;
finalctag=reshape(ctag(:,20),dimval,dimval);
initialctag=ctag(:,1);
finalgridconc=statevars{20}(1:3*dimval*dimval);
stateV=statevars{20}(3*dimval*dimval+1:end);
stateV2=reshape(stateV,size(stateV,1)/8,8);
figure
imagesc(finalctag);
%% Set up Initial conditions
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

intrazeb1=1.125; % baseline was 1 but 1.125 leads to prolonged P state
intrasnail2=1;
kdae_scale=1;
comment=['prolifer8'];


% Initial Seeding Pattern

% Use ctag from a simulation ****************
Fulltag=finalctag;
Fulltag=reshape(Fulltag,dimval*dimval,1);
initmatrix=Fulltag;


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
%cellvars=stateV2;
Fullstate=stateV2;
cellvars=Fullstate;
cellvars=reshape(cellvars,max(max(initmatrix))*size(icstate,1),1);

% Initial seeding matrix concentrations
imatstate=[0.1571,0.0072,0.01,2.4569*10^-5];

% parammatrix
parammatrix=modelparameters(:,1:max(max(finalctag)));

griddim=dimval*dimval;
% Alternatively, get grid concentrations from simulation, (Could use
% epithelial or mesenchymal grid here) **********************
matrvars=nan(dimval*dimval,3);
matrvars(:,1)=finalgridconc(1:1*dimval*dimval);
matrvars(:,2)=finalgridconc(1*dimval*dimval+1:2*dimval*dimval);
matrvars(:,3)=finalgridconc(2*dimval*dimval+1:3*dimval*dimval);
matrvars=reshape(matrvars,3*dimval*dimval,1);


%%
clc

run_single_cpm_fem_ode(dimval,mcsval,tgfbval,spreadval,tmaxval,Jhalf,dividing,extrascale1,extrascale2,intrazeb1, intrasnail2,pZEB,pR200,kdae_scale,comment,parammatrix,initmatrix,cellvars,matrvars);
pause(1)
close all

