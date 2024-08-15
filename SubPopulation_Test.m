Emat='test_100p_3000mcs_0tgfb_0.06spread1.25div_extraA0.96_extraBT0.6_ZEB1.125_SNAIL1_kDAE.scale1_Egrid.mat';
Mmat='test_100p_2500mcs_1tgfb_0.06spread1.25div_extraA0.96_extraBT0.6_ZEB1.125_SNAIL1_kDAE.scale1_Egrid.mat';
load('Cell_types.mat');
%% Epithelial subpopulation test- EMT Prone
load(Emat) %_inelasticity_1000
dimval=100;
finalctag=reshape(ctag(:,20),dimval,dimval);
initialctag=ctag(:,1);
finalgridconc=statevars{20}(1:3*dimval*dimval);
stateV=statevars{20}(3*dimval*dimval+1:end);
stateV2=reshape(stateV,size(stateV,1)/8,8);
figure
imagesc(finalctag);
%%
%ECMbonusrange=[1/1.1 1/1.25 1/1.5];
ECMbonusrange=[1/(4/3) 0.7];
parfor ECMb=[1:2]
    ECMbonus=ECMbonusrange(ECMb);
mcsval=300;
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
extrascale1=0.96*ECMbonus; % ka- base is 0.96
extrascale2=0.6*ECMbonus; % ECM TGFB Binding Rate kBT- base is 0.6
intrazeb1=1.125; % baseline was 1 but 1.125 leads to prolonged P state
intrasnail2=1; % base is 1
kdae_scale=1/ECMbonus;  % base is 1
ke_scale=1*ECMbonus;
%P_count=pval;

comment=['EMTprone_ECMscale' num2str(ECMbonus)];


% Initial Seeding Pattern

% Use ctag from a simulation ****************
Fulltag=finalctag;
Fulltag=reshape(Fulltag,dimval*dimval,1);
initmatrix=Fulltag;


%Initial Seeding cell concentrations
icstate = [.01; .01; .38; ...
    .03; .01; .35; 3.2; 0];  % initial state, endoTGFB 0.16 omitted, NEED TO KEEP THIS

% Use cell final statevars from simulation (Could use epithelial or
% mesenchymal here)
%cellvars=stateV2;
Fullstate=stateV2;
cellvars=Fullstate;
cellvars=reshape(cellvars,max(max(initmatrix))*size(icstate,1),1);

% Initial seeding matrix concentrations
imatstate=[0.1571,0.0072,0.01,2.4569*10^-5];

% parammatrix
dat1=datasample([1:1:100],max(max(finalctag)),'Replace',true);
parammatrix=P_EMT_prone(dat1,:)'; %modelparameters(:,1:max(max(finalctag)));
%
griddim=dimval*dimval;
% Alternatively, get grid concentrations from simulation, (Could use
% epithelial or mesenchymal grid here) **********************
matrvars=nan(dimval*dimval,3);
matrvars(:,1)=finalgridconc(1:1*dimval*dimval);
matrvars(:,2)=finalgridconc(1*dimval*dimval+1:2*dimval*dimval);
matrvars(:,3)=finalgridconc(2*dimval*dimval+1:3*dimval*dimval);
matrvars=reshape(matrvars,3*dimval*dimval,1);
%
clc
run_single_cpm_fem_ode(dimval,mcsval,tgfbval,spreadval,tmaxval,Jhalf,dividing,extrascale1,extrascale2,intrazeb1, intrasnail2,pZEB,pR200,kdae_scale,ke_scale,comment,parammatrix,initmatrix,cellvars,matrvars);
pause(1)
close all
end
%% Epithelial subpopulation test- EMT Resistant
load(Emat) %_inelasticity_1000
dimval=100;
finalctag=reshape(ctag(:,20),dimval,dimval);
initialctag=ctag(:,1);
finalgridconc=statevars{20}(1:3*dimval*dimval);
stateV=statevars{20}(3*dimval*dimval+1:end);
stateV2=reshape(stateV,size(stateV,1)/8,8);
figure
imagesc(finalctag);
%%
ECMbonusrange=[1.155 1.165 1.17];
load('test_50p_2500mcs_0tgfb_0.075spread1.25div_extraA0.96_extraBT0.6_ZEB1.125_SNAIL1_kDAE.scale1_base.mat');
matrvars=nan(dimval*dimval,3);
finalgridconc50=statevars{end}(1:3*50*50);
matrvars(:,1)=finalgridconc50(1+1225);
matrvars(:,2)=finalgridconc50(2501+1225);
matrvars(:,3)=finalgridconc50(5001+1225);
matrvars=reshape(matrvars,3*dimval*dimval,1);
for tgfbval=[0 0.5]
parfor ECMb=[1:3]
    ECMbonus=ECMbonusrange(ECMb);
mcsval=800;
%ECMbonus=1.25;
spreadval=0.015*4;
tmaxval=0;
Jmin=4819*15;
dividing=0; %1.25;%/4; PROLIFERATION RATE TURNED OFF????
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
extrascale1=0.96*ECMbonus; % ka- base is 0.96
extrascale2=0.6*ECMbonus; % ECM TGFB Binding Rate kBT- base is 0.6
intrazeb1=1.125; % baseline was 1 but 1.125 leads to prolonged P state
intrasnail2=1; % base is 1
kdae_scale=1/ECMbonus;  % base is 1
ke_scale=1*ECMbonus;
% pval=1;
% P_count=pval;


% Initial Seeding Pattern

% Use ctag from a simulation ****************
Fulltag=finalctag;
Fulltag=reshape(Fulltag,dimval*dimval,1);
initmatrix=Fulltag;


%Initial Seeding cell concentrations
icstate = [.01; .01; .38; ...
    .03; .01; .35; 3.2; 0];  % initial state, endoTGFB 0.16 omitted, NEED TO KEEP THIS

% Use cell final statevars from simulation (Could use epithelial or
% mesenchymal here)
%cellvars=stateV2;
Fullstate=stateV2;
cellvars=Fullstate;
cellvars=reshape(cellvars,max(max(initmatrix))*size(icstate,1),1);

% Initial seeding matrix concentrations
imatstate=[0.1571,0.0072,0.01,2.4569*10^-5];
% Pmin=(P_count-1)*max(max(finalctag))+1;
% Pmax=P_count*max(max(finalctag));
comment=['EMTres' '_ECMscale' num2str(ECMbonus)];
% parammatrix
dat1=datasample([1:1:95],max(max(finalctag)),'Replace',true);
parammatrix=P_EMT_res(dat1,:)'; %modelparameters(:,1:max(max(finalctag)));

griddim=dimval*dimval;
% Alternatively, get grid concentrations from simulation, (Could use
% epithelial or mesenchymal grid here) **********************
% matrvars=nan(dimval*dimval,3);
% matrvars(:,1)=finalgridconc(1:1*dimval*dimval);
% matrvars(:,2)=finalgridconc(1*dimval*dimval+1:2*dimval*dimval);
% matrvars(:,3)=finalgridconc(2*dimval*dimval+1:3*dimval*dimval);
% matrvars=reshape(matrvars,3*dimval*dimval,1);
%
clc
disp(comment)
run_single_cpm_fem_ode(dimval,mcsval,tgfbval,spreadval,tmaxval,Jhalf,dividing,extrascale1,extrascale2,intrazeb1, intrasnail2,pZEB,pR200,kdae_scale,ke_scale,comment,parammatrix,initmatrix,cellvars,matrvars);
pause(1)
close all
end
end
%% Epithelial subpopulation test- EMT Partial-Mesenchymal Resistant
load(Emat) %_inelasticity_1000
dimval=100;
finalctag=reshape(ctag(:,20),dimval,dimval);
initialctag=ctag(:,1);
finalgridconc=statevars{20}(1:3*dimval*dimval);
stateV=statevars{20}(3*dimval*dimval+1:end);
stateV2=reshape(stateV,size(stateV,1)/8,8);
figure
imagesc(finalctag);
%%
for ECMbonus=[5]
    parfor pval=1:6
mcsval=900;
tgfbval=1;

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
extrascale1=0.96*ECMbonus; % ka- base is 0.96
extrascale2=0.6*ECMbonus; % ECM TGFB Binding Rate kBT- base is 0.6
intrazeb1=1.125; % baseline was 1 but 1.125 leads to prolonged P state
intrasnail2=1; % base is 1
kdae_scale=1/ECMbonus;  % base is 1
ke_scale=1*ECMbonus;
P_count=pval;




% Initial Seeding Pattern

% Use ctag from a simulation ****************
Fulltag=finalctag;
Fulltag=reshape(Fulltag,dimval*dimval,1);
initmatrix=Fulltag;


%Initial Seeding cell concentrations
icstate = [.01; .01; .38; ...
    .03; .01; .35; 3.2; 0];  % initial state, endoTGFB 0.16 omitted, NEED TO KEEP THIS

% Use cell final statevars from simulation (Could use epithelial or
% mesenchymal here)
%cellvars=stateV2;
Fullstate=stateV2;
cellvars=Fullstate;
cellvars=reshape(cellvars,max(max(initmatrix))*size(icstate,1),1);

% Initial seeding matrix concentrations
imatstate=[0.1571,0.0072,0.01,2.4569*10^-5];

% parammatrix
Pmin=(P_count-1)*max(max(finalctag))+1;
Pmax=P_count*max(max(finalctag));
comment=['EMT-PMres' num2str(Pmin) 'to' num2str(Pmax) 'ECMscale' num2str(ECMbonus)];
parammatrix=P_PM_res(Pmin:Pmax,:)'; %modelparameters(:,1:max(max(finalctag)));

griddim=dimval*dimval;
% Alternatively, get grid concentrations from simulation, (Could use
% epithelial or mesenchymal grid here) **********************
matrvars=nan(dimval*dimval,3);
matrvars(:,1)=finalgridconc(1:1*dimval*dimval);
matrvars(:,2)=finalgridconc(1*dimval*dimval+1:2*dimval*dimval);
matrvars(:,3)=finalgridconc(2*dimval*dimval+1:3*dimval*dimval);
matrvars=reshape(matrvars,3*dimval*dimval,1);
%
clc
run_single_cpm_fem_ode(dimval,mcsval,tgfbval,spreadval,tmaxval,Jhalf,dividing,extrascale1,extrascale2,intrazeb1, intrasnail2,pZEB,pR200,kdae_scale,ke_scale,comment,parammatrix,initmatrix,cellvars,matrvars);
pause(1)
close all

    end
end

%% Mesenchymal subpopulation test
% Load Simulation- Mesenchymal
% get ctag matrix and intra/extracellular concentrations
load(Mmat) %_inelasticity_1000
dimval=100;
finalctag=reshape(ctag(:,end),dimval,dimval);
initialctag=ctag(:,1);
finalgridconc=statevars{end}(1:3*dimval*dimval);


halfc=10/2; % Cluster size is 10+1, must be an integer


% finalctag=reshape(finalctag,100,100);
% Visualize Initial Grid
loc=finalctag; % Rotates matrix to match video orientation
figure
imagesc(loc)
cells=max(max(loc));
%Centroid stuff
for i=1:cells
    [f1,f2]=find(loc == i);
    hold on
    xc=round(mean(f2),0);
    yc=round(mean(f1),0);
    text(xc,yc,[num2str(i)],'Color','w')
end

% Select cells
% [y,x] = ginput;
% %
% selected=zeros(size(x));
% for i=1:size(x,1)
%     x1=round(x(i),0);
%     y1=round(y(i),0);
% selected(i)=loc(x1,y1);
% end
halfval=dimval/2;
cluster1=loc(halfval-halfc:halfval+halfc,halfval-halfc:halfval+halfc);
unq2=unique(cluster1);
selected=unq2;
% Remove everything but select cells
remain=selected;
finaltag2=finalctag;
for pix=1:dimval*dimval
    check1=0;
    for r=1:size(remain,1)
        if(finaltag2(pix)==remain(r))
            check1=1;
        end
    end
    if check1==0
        finaltag2(pix)=0;
    end
    
end
% figure
% imagesc(finaltag2);
finaltag3=finaltag2;
% reorder ctag to remove gaps in cell # order
changed=0;
finaltag4=zeros(dimval,dimval);
testf1=finaltag3==remain(r);
testf2=zeros(dimval,dimval);
testf2(finaltag3==remain(r))=100;

    for r=1:size(remain,1)  
        finaltag4(finaltag3==remain(r))=r;
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

figure
imagesc(finaltag4)

%Centroid stuff
for i=1:NRc
    [f1,f2]=find(finaltag4 == i);
    hold on
    xc=round(mean(f2),0);
    yc=round(mean(f1),0);
    disp(['Marked: ' num2str(i)]);
    text(xc,yc,[num2str(i)],'Color','w')
end
Mtag=finaltag4;
% Give state variables to each cell
stateV=statevars{end}(3*dimval*dimval+1:end);
stateV2=reshape(stateV,size(stateV,1)/8,8);
Mstate=stateV2(remain,:);
%%
load('test_50p_2500mcs_0tgfb_0.075spread1.25div_extraA0.96_extraBT0.6_ZEB1.125_SNAIL1_kDAE.scale1_base.mat');
matrvars=nan(dimval*dimval,3);
finalgridconc50=statevars{end}(1:3*50*50);
matrvars(:,1)=finalgridconc50(1+1225);
matrvars(:,2)=finalgridconc50(2501+1225);
matrvars(:,3)=finalgridconc50(5001+1225);
matrvars=reshape(matrvars,3*dimval*dimval,1);

for ECMbonus=[1/1.8 1/1.75]
parfor pval=1:6
P_count=pval+10;
mcsval=600;
tgfbval=1;

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
extrascale1=0.96*ECMbonus; % ka- base is 0.96
extrascale2=0.6*ECMbonus; % ECM TGFB Binding Rate kBT- base is 0.6
intrazeb1=1.125; % baseline was 1 but 1.125 leads to prolonged P state
intrasnail2=1; % base is 1
kdae_scale=1/ECMbonus;  % base is 1
ke_scale=1*ECMbonus;


% Initial Seeding Pattern

% Use ctag from a simulation ****************
Fulltag=Mtag;
Fulltag=reshape(Fulltag,dimval*dimval,1);
initmatrix=Fulltag;


%Initial Seeding cell concentrations
icstate = [.01; .01; .38; ...
    .03; .01; .35; 3.2; 0];  % initial state, endoTGFB 0.16 omitted, NEED TO KEEP THIS

% Use cell final statevars from simulation (Could use epithelial or
% mesenchymal here)
%cellvars=stateV2;
Fullstate=Mstate;
cellvars=Fullstate;
cellvars=reshape(cellvars,max(max(initmatrix))*size(icstate,1),1);

% Initial seeding matrix concentrations
imatstate=[0.1571,0.0072,0.01,2.4569*10^-5];

Pmin=(P_count-1)*max(max(Fulltag))+1;
Pmax=P_count*max(max(Fulltag));
comment=['Egrid_METprone' num2str(Pmin) 'to' num2str(Pmax) 'ECMscale' num2str(ECMbonus)];
% parammatrix
parammatrix=P_MET_prone(Pmin:Pmax,:)'; %modelparameters(:,1:max(max(finalctag)));

griddim=dimval*dimval;
% Alternatively, get grid concentrations from simulation, (Could use
% epithelial or mesenchymal grid here) **********************
% finalgridconc=finalgridconcE;
% matrvars=nan(dimval*dimval,3);
% matrvars(:,1)=finalgridconc(1:1*dimval*dimval);
% matrvars(:,2)=finalgridconc(1*dimval*dimval+1:2*dimval*dimval);
% matrvars(:,3)=finalgridconc(2*dimval*dimval+1:3*dimval*dimval);
% matrvars=reshape(matrvars,3*dimval*dimval,1);

%
clc
run_single_cpm_fem_ode(dimval,mcsval,tgfbval,spreadval,tmaxval,Jhalf,dividing,extrascale1,extrascale2,intrazeb1, intrasnail2,pZEB,pR200,kdae_scale,ke_scale,comment,parammatrix,initmatrix,cellvars,matrvars);
pause(1)
close all
end
end

%%
run_single_cpm_fem_ode(dimval,mcsval,tgfbval,spreadval,tmaxval,Jhalf,dividing,extrascale1,extrascale2,intrazeb1, intrasnail2,pZEB,pR200,kdae_scale,ke_scale,comment,parammatrix,initmatrix,cellvars,matrvars);
