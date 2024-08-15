%% PMN Scenarios 
% base grids for creating mask
Emat='test_100p_2501mcs_0tgfb_0.075spread1.25div_extraA0.96_extraBT0.6_ZEB1.125_SNAIL1_kDAE.scale1_base.mat';
Mmat='test_100p_2501mcs_1tgfb_0.075spread1.25div_extraA0.96_extraBT0.6_ZEB1.125_SNAIL1_kDAE.scale1_base.mat';
gridmat='test_50p_2500mcs_0tgfb_0.075spread1.25div_extraA0.96_extraBT0.6_ZEB1.125_SNAIL1_kDAE.scale1_base.mat';
load('Cell_types.mat');
load('baseline_parameters');
%% Create epithelial periphery
load(Emat) 
dimval=100;
finalctag=reshape(ctag(:,end),dimval,dimval);

% Visualize Initial Grid
loc=finalctag; % Rotates matrix to match cell location to video orientation
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

%% Create mesenchymal center- 1st half
% get ctag matrix and intra/extracellular concentrations
load(Emat) %_inelasticity_1000
dimval=100;
finalctag=reshape(ctag(:,end),dimval,dimval);
initialctag=ctag(:,1);
%
halfc=8/2; % Cluster size is 10+1, must be an integer, so numerator must be even number


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


halfval=dimval/2;
%cluster1=loc(halfval-halfc:halfval+halfc,halfval-halfc:halfval+halfc);
cluster1=loc(halfval-halfc:halfval+halfc-1,halfval-halfc+1:halfval+halfc-1); %emat
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

%% Create mesenchymal center- 2nd half
% get ctag matrix and intra/extracellular concentrations
load(Mmat) %_inelasticity_1000
stateV=statevars{end}(3*dimval*dimval+1:end);
stateV2=reshape(stateV,size(stateV,1)/8,8);
halfremain=size(remain,1)/2;
Mstate(halfremain+1:end,:)=stateV2(1:halfremain,:);

%% Create ECM grid
load(gridmat);
%

grid=50;
fulltime=2500;
solTvar=zeros(grid*grid,1);
solEvar=zeros(grid*grid,1);
asmEvar=zeros(grid*grid,1);
ecmTvar=zeros(grid*grid,1);
for MCS=fulltime
maskG=ctag(:,MCS);
stateV=statevars{MCS};
NRC=size(csize{MCS},1);
k0T_val=zeros(1,NRC);
k0T_val(:)=0.06;
kT_val=zeros(1,NRC);
kT_val(:)=1.2;
nx=grid*grid;
MR200_val=zeros(sqrt(nx),sqrt(nx));
R2Val=stateV((grid*grid*3)+(5*NRC)+1: (grid*grid*3)+(6*NRC));
for i=1:NRC
    M200_temp=R2Val(i);
    MR200_val(find(maskG==i))=M200_temp;
end
MR200_val=reshape(MR200_val,nx,1);

JT_val=zeros(1,NRC);
JT_val(:)=0.06;
nr2_val=zeros(1,NRC);
nr2_val(:)=2;
kT200 = k0T_val + kT_val./(1 + ([MR200_val]./JT_val).^nr2_val);
kT200(isnan(kT200))=0;

solEv=stateV(1:(grid*grid));
totalEv=stateV((grid*grid+1):(grid*grid*2));
totalTv=stateV((grid*grid*2+1):(grid*grid*3));
kbT=10.656*0.6;

Ncad_val=zeros(sqrt(nx),sqrt(nx));
NVal=stateV((grid*grid*3)+(7*NRC)+1: (grid*grid*3)+(8*NRC));
for i=1:NRC
    Ncad_temp=NVal(i);
    Ncad_val(find(maskG==i))=Ncad_temp;
end
Ncad_val=reshape(Ncad_val,nx,1);

param_kdae=1.225;
fdAE=5; % Scaling Factor
kdAE2=(720/500/25)*param_kdae;
kdAE1=kdAE2/fdAE; % from ECM degradation graph
kdAE=(kdAE2-kdAE1).*Ncad_val/3.1515 + kdAE1; % Uses Ncad so calculated here
kdAE=kdAE*1;

ecmTv=(kbT.*(totalTv+totalEv)+kdAE - sqrt( (kbT.*(totalTv+totalEv)+kdAE).^2 - 4.*kbT.^2.*totalTv.*totalEv))./(2*kbT);
solTv=totalTv-ecmTv;
asmEv=totalEv-ecmTv;

% storing variables
solTvar(:,1)=solTv;
solEvar(:,1)=solEv;
asmEvar(:,1)=asmEv;
ecmTvar(:,1)=ecmTv;

end

%
cellloc=maskG>0;
totalvals=[solTvar(:,1),solEvar(:,1),asmEvar(:,1),ecmTvar(:,1)];
midvals=totalvals(1225,:);
avgvals=mean(totalvals(cellloc,:))

finalgridconc=stateV(1:3*grid*grid);
matrvars=nan(dimval*dimval,3);
matrvars(:,1)=finalgridconc(1+1225);
matrvars(:,2)=finalgridconc(2501+1225);
matrvars(:,3)=finalgridconc(5001+1225);
matrvars=reshape(matrvars,3*dimval*dimval,1);
disp([num2str(finalgridconc(1+1225)) ', ' num2str(finalgridconc(2501+1225)) ', ' num2str(finalgridconc(5001+1225))]);

%% Combine everything to create mask
ENRc=max(max(Etag));
Mtag2=Mtag;

Mtag2(Mtag>0)=Mtag2(Mtag>0)+ENRc; % change Mtag numbers to be after Etag

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

Fullstate=[Estate; Mstate]; % state variables for each cell

%% Other parameters

tgfbval=0; % ExoTGFB added to the system

spreadval=0.015*4;
%tmaxval=0;
Jmin=4819*15; % 15 is base
dividing=0; %1.25/1000;%/4; % Proliferation rate (1.25 is base, but /50 makes proliferation very unlikely)
al=0.05;
n=2;
Jhalf0=Jmin/((al^-1-1)^(1/n));
Tmax=2*.15; % .15 is base
tmax0=1*.15;
Jmin= Jhalf0*(al^-1 -1) ^(1/n) / (2*Tmax/tmax0 - 1)^(1/n);
Jhalf=Jmin/((al^-1-1)^(1/n));

pZEB=0.8;
pR200=1.2;
tmaxval=Tmax;
intrazeb1=1.125; % baseline was 1 but 1.125 leads to prolonged P state
intrasnail2=1; % base is 1


% Initial Seeding Pattern

% Use ctag from initially created mask **************** 
Fulltag=reshape(Fulltag,dimval*dimval,1); % Fulltag is combined grid
initmatrix=Fulltag;


%Initial Seeding cell concentrations
icstate = [.01; .01; .38; ...
    .03; .01; .35; 3.2; 0];  % initial state, endoTGFB 0.16 omitted, NEED TO KEEP THIS

% Use cell final statevars from simulation (Could use epithelial or
% mesenchymal here)
%cellvars=stateV2;
cellvars=Fullstate;
cellvars=reshape(cellvars,max(max(initmatrix))*size(icstate,1),1);

% Initial seeding matrix concentrations
imatstate=[0.1571,0.0072,0.01,2.4569*10^-5];


griddim=dimval*dimval;
%% ECM parameters
parfor mc=1:5
    mcsval=3000+mc;
ECMrange=[1];

    % Parameters for cells
Esize=size(Estate,1);
Msize=size(Mstate,1);
epiparams=zeros(Esize,43);
mesparams=zeros(Msize,43);
Esample=datasample([1:1:size(P_EMT_res,1)],Esize,'Replace',true);
Msample=datasample([1:1:size(P_EMT_prone,1)],Msize,'Replace',true);
Msample2=datasample([1:1:size(P_MET_prone,1)],Msize,'Replace',true);
for e=1:Esize
epiparams(e,:)=P_EMT_res(Esample(e),:); % baseparameters ; %baseparameters ******
end
for m=1:Msize/2
mesparams(m,:)=P_EMT_prone(Msample(m),:);%P_MET_prone(Msample(m),:); % baseparameters P
end
for m=Msize/2+1:Msize
mesparams(m,:)=P_MET_prone(Msample2(m),:); % baseparameters P
end
parammatrix=[epiparams; mesparams]'; 

    
    ECMbonus= 1.1;%ECMrange(ec);% baseline scaling is 1. >1 is pro-EMT, <1 is anti-EMT

extrascale1=0.96*ECMbonus; % ka- base is 0.96
extrascale2=0.6*ECMbonus; % ECM TGFB Binding Rate kBT- base is 0.6
kdae_scale=1/ECMbonus;  % base is 1
ke_scale=1*ECMbonus;
comment=['ECMscale' num2str(ECMbonus) '__26EMTresPer_METprone&EMTproneCent_ConvertGrid'];

%
clc
run_single_cpm_fem_ode(dimval,mcsval,tgfbval,spreadval,tmaxval,Jhalf,dividing,extrascale1,extrascale2,intrazeb1, intrasnail2,pZEB,pR200,kdae_scale,ke_scale,comment,parammatrix,initmatrix,cellvars,matrvars);
pause(1)
end