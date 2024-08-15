Emat='test_100p_1501mcs_1tgfb_0.06spread_ECMscale1_PMN_TransformableBefore_0.015789_divf_0.1.mat';
Mmat='test_100p_2501mcs_1tgfb_0.075spread1.25div_extraA0.96_extraBT0.6_ZEB1.125_SNAIL1_kDAE.scale1_base.mat';
gridmat=Emat;
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



Etag=finalctag;
% Give state variables to each cell
stateV=statevars{end}(3*dimval*dimval+1:end);
stateV2=reshape(stateV,size(stateV,1)/8,8);
Estate=stateV2;

%% Create mesenchymal center
% get ctag matrix and intra/extracellular concentrations
%load(Emat) %_inelasticity_1000
dimval=100;
%finalctag=reshape(ctag(:,end),dimval,dimval);
initialctag=ctag(:,1);
blankgrid=zeros(100,100);
count=5;
for i=45:55
    for j=45:55
        if (count>0 && finalctag(i,j)==0)
            blankgrid(i,j)=max(max(finalctag))+count;
            count=count-1;
            disp('Tumor cell added')
        end
    end
end
%
figure
subplot(1,2,1)
imagesc(blankgrid)
subplot(1,2,2)
combined1=blankgrid+finalctag;
imagesc(combined1);


%% Create ECM grid
load(Emat);
%
grid=100;
fulltime=3000;
stateV=statevars{end};
finalgridconc=stateV(1:3*grid*grid);
matrvars=nan(dimval*dimval,3);
matrvars(:,1)=finalgridconc(1:10000);
matrvars(:,2)=finalgridconc(10001:20000);
matrvars(:,3)=finalgridconc(20001:30000);
matrvars=reshape(matrvars,3*dimval*dimval,1);

%% Combine everything to create mask
%ENRc=max(max(Etag));
load(Mmat);
stateV=statevars{end}(3*dimval*dimval+1:end);
stateV2=reshape(stateV,size(stateV,1)/8,8);
Mstate=stateV2(1:5,:);

Fulltag=combined1;

Fullstate=[Estate; Mstate]; % state variables for each cell


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
% Use ctag from initially created mask **************** 
Fulltag=reshape(Fulltag,dimval*dimval,1); % Fulltag is combined grid
initmatrix=Fulltag;

% Use cell final statevars from simulation (Could use epithelial or
% mesenchymal here)
%cellvars=stateV2;
cellvars=Fullstate;
cellvars=reshape(cellvars,max(max(initmatrix))*size(icstate,1),1);

% Initial seeding matrix concentrations
imatstate=[0.1571,0.0072,0.01,2.4569*10^-5];

griddim=dimval*dimval;

%% Cells pop_index
load(Emat);
%%
total1=NRc(end);
og=size(pop_indexOG,1);
adding=total1-og;
pop_indexNew=zeros(total1+5,1);
pop_indexNew(1:og)=pop_indexOG;
pop_indexNew(og+1:og+adding)=2;
pop_indexNew(og+adding+1:end)=3;

%% ECM parameters
ECMrange=[1 1/1.75 1.1];

ksval=1.9;
P_nontransform=P;
P_nontransform(7)=P_nontransform(7)/ksval;
ec=1;
    ECMbonus= ECMrange(ec);% baseline scaling is 1. >1 is pro-EMT, <1 is anti-EMT
    if ec==1
     gridname='Base';
%     elseif ec==2
%       gridname='Revert';
%     elseif ec==3
%         gridname='Convert';
     end
divfactor=0.1;
comment=['ECMscale' num2str(ECMbonus) '_PMN_TransformableAfter_' num2str(P_nontransform(7)) '_divf_' num2str(divfactor)];
extrascale1=0.96*ECMbonus; % ka- base is 0.96
extrascale2=0.6*ECMbonus; % ECM TGFB Binding Rate kBT- base is 0.6
kdae_scale=1/ECMbonus;  % base is 1
ke_scale=1*ECMbonus;

tgfrange=[0 0.1 1 5 10];%[0 0.1 0.2 0.3 0.4 0.5]; %min 0.1 needed for differentiation

prob_pop = [0.5 0.5];
NRcsimulation=total1+5;
%pop_index = zeros(NRcsimulation,1);
tempparam=zeros(NRcsimulation,43);
prolif=zeros(NRcsimulation,1);
for i1=1:size(pop_indexNew,1)
%pop_index(i1)=find(rand<cumsum(prob_pop),1,'first');
if(pop_indexNew(i1)==1)
tempparam(i1,:)=P_nontransform;
prolif(i1)=0;
elseif(pop_indexNew(i1)==2)
    tempparam(i1,:)=P;
    prolif(i1)=1.25*divfactor;
elseif(pop_indexNew(i1)==3)
    tempparam(i1,:)=P;
    prolif(i1)=1.25;
end
end

dividing=prolif;


%
for i=1
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

mcsval=500;
clc
run_single_cpm_fem_ode(dimval,mcsval,tgfbval,spreadval,tmaxval,Jhalf,dividing,extrascale1,extrascale2,intrazeb1, intrasnail2,pZEB,pR200,kdae_scale,ke_scale,comment,parammatrix,initmatrix,cellvars,matrvars);
pause(1)

end





