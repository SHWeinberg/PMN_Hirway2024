%% PMN Scenario with simulation of kupffer and stellate cells and then tumor cells added
% base grids for creating mask
Emat='test_100p_2501mcs_0tgfb_0.075spread1.25div_extraA0.96_extraBT0.6_ZEB1.125_SNAIL1_kDAE.scale1_base.mat';
Mmat='test_100p_2501mcs_1tgfb_0.075spread1.25div_extraA0.96_extraBT0.6_ZEB1.125_SNAIL1_kDAE.scale1_base.mat';
%gridmat='test_50p_2500mcs_0tgfb_0.075spread1.25div_extraA0.96_extraBT0.6_ZEB1.125_SNAIL1_kDAE.scale1_base.mat';
gridmat='test_100p_1007mcs_0tgfb__ConfluenceK_divf_0.5_param_14_scale_0.5ECMscaled1_K.mat'; % Using confluent K non transformed cells as base

load('Cell_types.mat');
load('baseline_parameters');
%% Create tumor variables
load(Mmat);
stateV=statevars{end}(3*dimval*dimval+1:end);
stateV2=reshape(stateV,size(stateV,1)/8,8);
tumor=stateV2(1:5,:);

% Randomly seed cells as pixels on grid
dimval=100;
spreadval=0.015*4;
VOXSIZE = 2.5E-6; % [m]
CELLDIAM = 2.0E-5; % cell diameter [m]
CELLRVOX = CELLDIAM/2/VOXSIZE; % cell radius [pixels]
TARGETVOLUME = 3.1415 * CELLRVOX * CELLRVOX; % targetvolume [pixels]


randomtag=zeros(dimval,dimval);
spread=spreadval*2.5;

NRc = 0;
for vy = 1:dimval
    for vx = 1:dimval
        v = vx + (vy - 1) * dimval;
        if ((vx > 1) && (vx < dimval) && (vy > 1) && (vy < dimval)) % exclude outer rim
            if (rand < (spread/TARGETVOLUME)) %smaller value makes it more sparse- SUH 112619 OG=0.25
                NRc = NRc + 1;
                randomtag(v) = NRc;
                disp(['Added cell ' num2str(NRc)]);
                %                 ecadr(v) = ECAD_DENSITY;
            end
        end
    end
end
split=round(NRc/2,0)+1;
NRcsimulation=NRc;
imagesc(randomtag)
% Give state variables to each cell
load(Emat)
Ematgrid=sqrt(size(ctag,1));
stateV=statevars{end}(3*Ematgrid*Ematgrid+1:end);
stateV=reshape(stateV,size(stateV,1)/8,8);
avgstateV=mean(stateV);
%stateV=statevars{end}(3*dimval*dimval+1:end);
%stateV2=reshape(stateV,size(stateV,1)/8,8);
%Estate=stateV2(remain,:);
icstate = [.01; .01; .38; .03; .01; .35; 3.2; 0];  % initial state, endoTGFB 0.16 omitted, NEED TO KEEP THIS
stateV2=zeros(NRcsimulation,8);
for i=1:NRcsimulation
stateV2(i,:)=avgstateV;
end

% Create ECM grid
load(gridmat);
%

grid=sqrt(size(ctag,1));
fulltime=size(ctag,2);
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
matrvars(:,1)=finalgridconc(1+grid*grid/2-grid/2);
matrvars(:,2)=finalgridconc(1*grid*grid+1+grid*grid/2-grid/2);
matrvars(:,3)=finalgridconc(2*grid*grid+1+grid*grid/2-grid/2);
matrvars=reshape(matrvars,3*dimval*dimval,1);
disp([num2str(finalgridconc(1+grid*grid/2-grid/2)) ', ' num2str(finalgridconc(1*grid*grid+1+grid*grid/2-grid/2)) ', ' num2str(finalgridconc(2*grid*grid+1+grid*grid/2-grid/2))]);


% Other parameters

tgfbval=0; % ExoTGFB added to the system


%tmaxval=0;
Jmin=4819*15; % 15 is base
dividing=1.25*0.25; %1.25/1000;%/4; % Proliferation rate (1.25 is base, but /50 makes proliferation very unlikely)
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
Fulltag=randomtag;
% Use ctag from initially created mask **************** 
Fulltag=reshape(Fulltag,dimval*dimval,1); % Fulltag is combined grid
initmatrix=Fulltag;

% Use cell final statevars from simulation (Could use epithelial or
% mesenchymal here)
cellvars=stateV2;
%cellvars=Fullstate;
cellvars=reshape(cellvars,max(max(initmatrix))*size(icstate,1),1);

% Initial seeding matrix concentrations
imatstate=[0.1571,0.0072,0.01,2.4569*10^-5];


griddim=dimval*dimval;
% Creating full paramter table with intracellular and extracellular terms
Pfull=P;
Pfull(44)=10.656; % 0.008- calculated from 12th-14th paper extra.kbT
Pfull(45)=0.006; % from 10-40nm ECM fibronectin prod extra.kE1
Pfull(46)=0.0273; %extra.kE2
%mdlparams.kE=((kE2-kE1).*Ncad_val/3.1515 + kE1)*init.param_ke;
Pfull(47)=0.6; %from Tian Model extra.kdE
Pfull(48)=10; % Scaling Factor extra.fa
Pfull(49)=1/(Pfull(48)*12); %extra.ka1
Pfull(50)=(1/12); %*extra.param_ka; %extra.ka2 scaled by param.ka in tian_baselineparams_v4
%mdlparams.ka=(ka2-ka1).*Ncad_val/3.1515 + ka1;
Pfull(51)=5; % Scaling Factor %extra.fdAE
param_kdae=1.225;
Pfull(52)=(720/500/25)*param_kdae; %*extra.param_kdae; 
Pfull(53)=Pfull(52)/Pfull(51); % from ECM degradation graph extra.kdAE1
%mdlparams.kdAE=(kdAE2-kdAE1).*Ncad_val/3.1515 + kdAE1;
Pkupffer=Pfull;

% cell contact tables
JcmTab=cell(3,1);
JccTab=cell(3,3);
JcaTab=cell(3,3);

% EMT based Tables
% JcaTab{1}=[0.5 0.3 0.1; 0.3 0.175 0.05; 0.1 0.05 0];
% JcaTab{2}=[0.5 0.3 0.1; 0.3 0.175 0.05; 0.1 0.05 0];
% JcaTab{3}=[0.5 0.3 0.1; 0.3 0.175 0.05; 0.1 0.05 0];
% JcaTab{4}=[0.5 0.3 0.1; 0.3 0.175 0.05; 0.1 0.05 0];
% JcaTab{5}=[0.5 0.3 0.1; 0.3 0.175 0.05; 0.1 0.05 0];
% JcaTab{6}=[0.5 0.3 0.1; 0.3 0.175 0.05; 0.1 0.05 0];
% JcaTab{7}=[0.5 0.3 0.1; 0.3 0.175 0.05; 0.1 0.05 0];
% JcaTab{8}=[0.5 0.3 0.1; 0.3 0.175 0.05; 0.1 0.05 0];
% JcaTab{9}=[0.5 0.3 0.1; 0.3 0.175 0.05; 0.1 0.05 0];
% 
% JccTab{1}=[1.7500    1.7500    1.7500 ;    1.7500    1.4250    1.4250 ;    1.7500    1.4250    1.1000];
% JccTab{2}=[1.7500    1.7500    1.7500 ;    1.7500    1.4250    1.4250 ;    1.7500    1.4250    1.1000];
% JccTab{3}=[1.7500    1.7500    1.7500 ;    1.7500    1.4250    1.4250 ;    1.7500    1.4250    1.1000];
% JccTab{4}=[1.7500    1.7500    1.7500 ;    1.7500    1.4250    1.4250 ;    1.7500    1.4250    1.1000];
% JccTab{5}=[1.7500    1.7500    1.7500 ;    1.7500    1.4250    1.4250 ;    1.7500    1.4250    1.1000];
% JccTab{6}=[1.7500    1.7500    1.7500 ;    1.7500    1.4250    1.4250 ;    1.7500    1.4250    1.1000];
% JccTab{7}=[1.7500    1.7500    1.7500 ;    1.7500    1.4250    1.4250 ;    1.7500    1.4250    1.1000];
% JccTab{8}=[1.7500    1.7500    1.7500 ;    1.7500    1.4250    1.4250 ;    1.7500    1.4250    1.1000];
% JccTab{9}=[1.7500    1.7500    1.7500 ;    1.7500    1.4250    1.4250 ;    1.7500    1.4250    1.1000];
% 
% JcmTab{1}=[2; 1.5; 1];
% JcmTab{2}=[2; 1.5; 1];
% JcmTab{3}=[2; 1.5; 1];


JcaTab{1,1}=[0.3 0.2 0.1; 0.2 0.15 0.05; 0.1 0.05 0]; %k-k
JcaTab{2,2}=[0.3 0.2 0.1; 0.2 0.15 0.05; 0.1 0.05 0]; %s-s
JcaTab{3,3}=[0.5 0.3 0.1; 0.3 0.175 0.05; 0.1 0.05 0]; %t-t
JcaTab{1,2}=[0.3 0.2 0.1; 0.2 0.15 0.05; 0.1 0.05 0]; %k-s
JcaTab{3,1}=[0.5 0.6 0.7; 0.3 0.45 0.6; 0.1 0.25 0.4]; %t-k
JcaTab{3,2}=[0.5 0.45 0.4; 0.3 0.25 0.2; 0.1 0.05 0]; %t-s
JcaTab{2,1}=JcaTab{1,2}'; %s-k
JcaTab{1,3}=JcaTab{3,1}'; %k-t
JcaTab{2,3}=JcaTab{3,2}'; %s-t

JccTab{1,1}=[1.5    1.5    1.5 ;    1.5    1.3    1.25 ;    1.5    1.25    1.1]; %k-k
JccTab{2,2}=[1.5    1.5    1.5 ;    1.5    1.3    1.25 ;    1.5    1.25    1.1]; %s-s
JccTab{3,3}=[1.75    1.75    1.75 ;    1.75    1.4250    1.4250 ;    1.75    1.4250    1.1]; %t-t
JccTab{1,2}=[1.5    1.5    1.5 ;    1.5    1.3    1.25 ;    1.5    1.25    1.1]; %k-s
JccTab{3,1}=[1.5    1.45    1.4 ;    1.4    1.35    1.3 ;    1.3    1.25    1.2]; %t-k
JccTab{3,2}=[1.3    1.275    1.25 ;    1.25    1.2    1.15 ;    1.2    1.15    1]; %t-s

JccTab{1,3}=JccTab{3,1}'; %k-t
JccTab{2,1}=JccTab{1,2}'; %s-k
JccTab{2,3}=JccTab{3,2}'; %s-t



JcmTab{1}=[1.5; 1.25; 1]; %k
JcmTab{2}=[1.5; 1.25; 1]; %s
JcmTab{3}=[2; 1.5; 1];% t

%% ECM parameters




ECMrange=[1 1/1.75 1.1];

ksval=1.9;
P_nontransform=P;
P_nontransform(7)=P_nontransform(7)/ksval;
ec=1;
    ECMbonus= ECMrange(ec);% baseline scaling is 1. >1 is pro-EMT, <1 is anti-EMT

     gridname='Base';

extrascale1=0.96*ECMbonus; % ka- base is 0.96
extrascale2=0.6*ECMbonus; % ECM TGFB Binding Rate kBT- base is 0.6
kdae_scale=1/ECMbonus;  % base is 1
ke_scale=1*ECMbonus;

tgfrange=[0 1 5];%[0 0.1 0.2 0.3 0.4 0.5]; %min 0.1 needed for differentiation
popcounter=zeros(6,2);
c=1;
tempscalerange=[1.32:0.01:1.37];
tempscaleindex=1:6;
    %tempscale=tempscalerange(tempscaleindex);
    tempscale=0; %ECM scaling parameter
    scaleval=0.5; %Kupffer cell scaling parameter for differentiation
    %[1:0.02:1.06]
    
    %tempscale=1;
    popnum=[1];
    
%for 
    
prob_pop = [0.5 0.5];

pop_index = zeros(NRcsimulation,1);
pop_index(:)=1;
pop_index(1:6)=2;
 pop_index(7:12)=1;
tempparam=zeros(NRcsimulation,53);
prolif=zeros(NRcsimulation,1);
diff1=scaleval;
paramvalue=14;

for i1=1:size(pop_index,1)
%pop_index(i1)=find(rand<cumsum(prob_pop),1,'first'); %Sets pobilation based
%on probability
if(pop_index(i1)==1) % Kupffer
    tempparam(i1,:)=Pkupffer;
    tempparam(i1,paramvalue)= tempparam(i1,paramvalue)*(1+diff1); % Decreased TGFB sensitivity
    tempparam(i1,52)=tempparam(i1,52)*tempscale; %kdae2- less 
     tempparam(i1,53)=tempparam(i1,52)/tempparam(i1,51); % kdae1
     tempparam(i1,45)=tempparam(i1,45)*tempscale; %ke1- less 
     tempparam(i1,46)=tempparam(i1,46)*tempscale; %ke2- less 
     tempparam(i1,49)=tempparam(i1,49)*tempscale; %ka1- less 
     tempparam(i1,50)=tempparam(i1,50)*tempscale; %ka2- less
     celltype='K';
elseif(pop_index(i1)==2) %Stellate
    tempparam(i1,:)=Pfull;
    tempparam(i1,paramvalue)= tempparam(i1,paramvalue);
    celltype='S';
elseif (pop_index(i1)==3) % Kupffer
    tempparam(i1,:)=Pkupffer;
    
    tempparam(i1,52)=tempparam(i1,52)*tempscale; %kdae2- less 
     tempparam(i1,53)=tempparam(i1,52)/tempparam(i1,51); % kdae1
     tempparam(i1,45)=tempparam(i1,45)*tempscale; %ke1- less 
     tempparam(i1,46)=tempparam(i1,46)*tempscale; %ke2- less 
     tempparam(i1,49)=tempparam(i1,49)*tempscale; %ka1- less 
     tempparam(i1,50)=tempparam(i1,50)*tempscale; %ka2- less
     celltype='T';
end
end
%popcounter(c,:)=[sum(pop_index(:)==1),sum(pop_index(:)==2)];
%c=c+1;
pause(1)
%dividing=prolif;
pop_table=struct;
%pop_table.div=[0 1.25];
divfactorS=0.5;
divfactorK=0.01;
%pop_table.div=[0 divfactor];
pop_table.div=[divfactorK divfactorS divfactorS];


% Other parameters
othermod=struct;
othermod.tmaxval=tmaxval;
othermod.Jhalf=Jhalf;
othermod.extrascale1=extrascale1;
othermod.extrascale2=extrascale2;
othermod.intrazeb1=intrazeb1;
othermod.intrasnail2=intrasnail2;
othermod.kdae_scale=kdae_scale;
othermod.ke_scale=ke_scale;
othermod.pZEB=pZEB;
othermod.pR200=pR200;
othermod.JcaTab=JcaTab;
othermod.JccTab=JccTab;
othermod.JcmTab=JcmTab;


%num2str(P_nontransform(7));
celltype='KST';
comment=['_BeforeTumor_divf_' num2str(divfactorS) '_param_' num2str(paramvalue) '_scale_' num2str(diff1) 'ECMscaled' num2str(tempscale) '_' num2str(IDval)];

%
i=3;
    % Parameters for cells
    tgfbval=tgfrange(i); % ExoTGFB added to the system
%Esize=size(NRc-split+1,1);
%Msize=size(NRc-split,1);
% lessparams=zeros(split-1,43);
% moreparams=zeros(NRcsimulation-split+1,43);
% Esample=datasample([1:1:size(P_EMT_res,1)],size(lessparams,1),'Replace',true);
% Msample=datasample([1:1:size(P_MET_prone,1)],size(moreparams,1),'Replace',true);
% for l=1:size(lessparams,1)
% lessparams(l,:)=P_nontransform;%P_EMT_res(Esample(e),:); % baseparameters ; %baseparameters ******
% end
% 
% for m=1:size(moreparams,1)
% moreparams(m,:)=P;%P_MET_prone(Msample(m),:); % baseparameters P
% end
parammatrix=[tempparam]'; 
%%
mcsval=15;
clc
run_single_cpm_fem_ode(dimval,mcsval,tgfbval,spreadval,comment,pop_index,pop_table,parammatrix,initmatrix,cellvars,matrvars,othermod);
pause(1)
