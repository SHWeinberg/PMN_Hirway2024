%% CPM Trial 5.1- Metastasis
% Have some epithelial cells with a cluster of mesenchymal cells invading
% the area

%% Load Simulation- Epithelial
% get ctag matrix and intra/extracellular concentrations
Emat='test_100p_2000mcs_0tgfb_0.3_tmax_0.06spread_Jhalf9574.3821_extraA0.96_extraBT0.6_ZEBscale0.8_R200scale1.2_ZEB1_SNAIL1_kDAE.scale1_size0.5.mat';
Mmat='test_100p_2000mcs_5tgfb_0.3_tmax_0.06spread_Jhalf9574.3821_extraA0.96_extraBT0.6_ZEBscale0.8_R200scale1.2_ZEB1_SNAIL1_kDAE.scale1_size0.5.mat';
load(Emat) %_inelasticity_1000
dimval=100;
finalctag=reshape(ctag(:,end),dimval,dimval);
initialctag=ctag(:,1);
finalgridconc=statevars{end}(1:3*dimval*dimval);
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

%Select Cells- ginput is by clicking, might take too much time
%[y,x] = ginput;
boundary=zeros(dimval-2,4);
boundary(:,1)=loc(2,2:dimval-1);
boundary(:,2)=loc(dimval-1,2:dimval-1);
boundary(:,3)=loc(2:dimval-1,2);
boundary(:,4)=loc(2:dimval-1,dimval-1);
boundary=reshape(boundary,(dimval-2)*4,1);
unq=unique(boundary);
%
% selected=zeros(size(x));
% for i=1:size(x,1)
%     x1=round(x(i),0);
%     y1=round(y(i),0);
% selected(i)=loc(x1,y1);
% end
% Remove everything but select cells
selected=unq;
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
Etag=finaltag4;
% Give state variables to each cell
stateV=statevars{end}(3*dimval*dimval+1:end);
stateV2=reshape(stateV,size(stateV,1)/8,8);
Estate=stateV2(remain,:);

%% Load Simulation- Mesenchymal
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


%% Combine grids
ENRc=max(max(Etag));
Mtag2=Mtag;

Mtag2(Mtag>0)=Mtag2(Mtag>0)+ENRc;

Fulltag=Etag+Mtag2;
figure
imagesc(Fulltag);
for i=1:max(max(Fulltag))
    [f1,f2]=find(Fulltag == i);
    hold on
    xc=round(mean(f2),0);
    yc=round(mean(f1),0);
    disp(['Marked: ' num2str(i)]);
    text(xc,yc,[num2str(i)],'Color','w')
end

Fullstate=[Estate; Mstate];

%%
load(Emat) %_inelasticity_1000
dimval=100;
finalgridconc=statevars{end}(1:3*dimval*dimval);
disp('E grid loaded');
%% Set up Initial conditions
mcsval=150;
tgfbval=0;
spreadval=0.015*4;
tmaxval=0;
Jmin=4819*15;
dividing=1.25/20;%/4;
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
comment=['PMN'];


% Initial Seeding Pattern

% Use ctag from a simulation ****************
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
cellvars=Fullstate;
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
% epithelial or mesenchymal grid here) **********************
matrvars=nan(dimval*dimval,3);
matrvars(:,1)=finalgridconc(1:1*dimval*dimval);
matrvars(:,2)=finalgridconc(1*dimval*dimval+1:2*dimval*dimval);
matrvars(:,3)=finalgridconc(2*dimval*dimval+1:3*dimval*dimval);
matrvars=reshape(matrvars,3*dimval*dimval,1);

% %% Change grid conc for scratch
% for i=1:10000
%     if(initmatrix(i)==0)
%         matrvars(i)=epigridconc(i);
%         matrvars(10000+i)=epigridconc(10000+i);
%         matrvars(20000+i)=epigridconc(0000+i);
%     end
% end
%%
clc

run_single_cpm_fem_ode(dimval,mcsval,tgfbval,spreadval,tmaxval,Jhalf,dividing,extrascale1,extrascale2,intrazeb1, intrasnail2,pZEB,pR200,kdae_scale,comment,initmatrix,cellvars,matrvars);
pause(1)
close all

