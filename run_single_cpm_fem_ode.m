function []=run_single_cpm_fem_ode(dimval,mcsval,tgfbval,spreadval,comment,pop_index,pop_table,parammatrix,initmatrix,cellvars,matrvars,othermod)
%v5.2
tic
% CPM-FEM v1 parameters
def.fmaflag=1;  % 0 - single, 1 - multicell FMA
def.ind=1;
% def.amp=2;
% def.ratio=2;
def.scurr=rng(othermod.rngseed);
% def.divthresh=0.005;
def.NVX = dimval;    % dimensions, pixels OG: 100x100
def.NRINC = mcsval;   % number of MCS steps OG: 100

def.NULL  = 0;
def.FALSE  = 0;
def.TRUE = 1;

def.NVY = def.NVX; 
def.NV = def.NVX * def.NVY;
def.NNX = def.NVX + 1; 
def.NNY = def.NNX; 
def.NN  = def.NNX * def.NNY;
def.NDOF = 2 * def.NN;
def.VOXSIZE = 2.5E-6; % [m]
def.MAXNRITER = 1000;
def.ACCURACY = .00001;

% material properties
def.YOUNGS = 10E3; % [Pa]
def.POISSON = .45; %

% loading
def.LOAD = 0; % 3E6
def.FORCE = (def.LOAD*def.VOXSIZE);

% cells
def.IMMOTILITY = 1.0; %50
def.CELLDIAM = 2.0E-5; % cell diameter [m]
def.CELLRVOX = def.CELLDIAM/2/def.VOXSIZE; % cell radius [pixels]
def.TARGETVOLUME = 3.1415 * def.CELLRVOX * def.CELLRVOX; % targetvolume [pixels]
def.INELASTICITY = 500.0*2; % [-] 1.0E20 % [/m4] % OG was 500.0

def.NOSTICKJ = 500000.0; %10000 [/m] contact penalty for non-adhesive surface
% def.JCM = def.amp * def.NOSTICKJ * def.VOXSIZE;  % cell-medium contact energy
% def.JCC = def.ratio * def.JCM; % cell-cell contact energy

def.MAXDHSTR = 10.0; % unscaled at the moment,
def.THRESHOLDSTIFF = 15E3;    % threshold stiffness for durotaxis
def.STIFFSENSITIVITY = .0005; % steepness of durotaxis sigmoid
def.STIFFENINGSTIFF = .1; % steepness of strain-stiffening

def.THICKNESS = 10E-6; % effective thickness of substrate
def.CELLFORCE = 1.0E-5/def.THICKNESS; % [N]

def.SQ05 = .707107; %sqrt(.5), used often enough to make this convenient
def.RAND_MAX = 32767;


% v2 parameters
% Initial state of cells defined here SUH 111220
cellmod.nstate = 8;  % number of state variables per cell % Change to 8 since endoTgfb is not used SUH 111320
cellmod.icstate = [.01; .01; .38; ...
    .03; .01; .35; 3.2; 0];  % initial state, endoTGFB 0.16 omitted
% Need to change IC state - SUH 111320
%Extracellular variable initial states
%cellmod.imatstate=[0.1571,0.01, 2.4569*10^-5,0.0072]; %solT solE asmE ecmT
cellmod.imatstate=[[0.140495366812545,0.00535985803143813,7.07650148176848e-05,0.00430341315227013]]; %solT solE asmE ecmT- based on epithelial monolayer from CPM3.1 - SUH 072121

% init.solT=0.1571;
% init.ecmT=0.0072;
% init.solE=0.01;
% init.asmE=2.4569*10^-5;
cellmod.initmatrix=initmatrix;
cellmod.cellvars=cellvars;
cellmod.matrvars=matrvars;

%**************************
% Population Based Cell Parameter Assignment
numcells=size(pop_index,1);
dividingL=zeros(numcells,1);
dividingH=zeros(numcells,1);
for i=1:numcells
    dividingL(i)=pop_table.divL(pop_index(i));
    dividingH(i)=pop_table.divH(pop_index(i));
end

%**************************

cellmod.params.method = 'MatODE'; % SimpleEuler, Euler, MatODE
cellmod.params.dt = 4.8;        % integration time step, min (for SimpleEuler, Euler methods)
cellmod.params.odefun = @CellFunc_v5;%@emt_tian_2013_v3; %Changed from OG script to v5- SUH 072721
cellmod.params.mcs_to_time = 4.8;  % min/MCS (for SimpleEuler, should be = to dt)

tmaxval=othermod.tmaxval;
Jhalf=othermod.Jhalf;
model_param=tian_model_baseline_params(tmaxval,Jhalf);
cellmod.tmaxval=othermod.tmaxval;
cellmod.Jhalf=othermod.Jhalf;
cellmod.extrascale1=othermod.extrascale1;
cellmod.extrascale2=othermod.extrascale2;
cellmod.intrazeb1=othermod.intrazeb1;
cellmod.intrasnail2=othermod.intrasnail2;
cellmod.kdae_scale=othermod.kdae_scale;
cellmod.ke_scale=othermod.ke_scale;
cellmod.pZEB=othermod.pZEB;
cellmod.pR200=othermod.pR200;
cellmod.JcaTab=othermod.JcaTab;
cellmod.JccTab=othermod.JccTab;
cellmod.JcmTab=othermod.JcmTab;

cellmod.parammatrix=parammatrix; %New parameter matrix from Mario's population study SUH 072221
% cell type based
cellmod.pop_table=pop_table;
cellmod.pop_index=pop_index;


%tian_model_baseline_params;  % OG script that defines model_param
%cellmod.params.model_params = model_param;
%cellmod.params.model_params= all_params;
%cellmod.params.model_params.scale = 3;   % scaling factor for speed up (>1) or slow down (< 1) model dynamics
%Scale changed from 7 to 3 - SUH 120519
% phenotype calculation parameters
cellmod.pheno.method = 'NcadConc';
cellmod.pheno.ncad_max = 3.1515;
cellmod.spread=spreadval; %seeding sparsity, lower value is more sparse
%cellmod.split=split;

% phenotype-dependent parameters 
scalef=0.5; % Changed from 2 to 0.5 SUH 7/12/21
phenoparams.targetvolE = 50.2640*scalef; %scaled up SUH 1/29/20
phenoparams.targetvolM = phenoparams.targetvolE*3;
phenoparams.pdivideE = 0.003*dividingL; % SUH 3/24/20- OG was 0.003
phenoparams.pdivideM = 0.003*dividingH; %phenoparams.pdivideE/3;
phenoparams.jcmE = 2.5;
phenoparams.jcmM = 2.5;
ratioE = 2; ratioM = 0.5;
phenoparams.jccE = ratioE*phenoparams.jcmE;
phenoparams.jccM = ratioM*phenoparams.jcmM;

% cell-dependent model parameters
cell_data.TGFb_max = 4;
cell_data.TGFBval=tgfbval; % exo-TGFB value here- SUH-112919

if 1
    cell_data.method = 'TGFb_Conn';
    cell_data.neighbor_max = 8;
elseif 0
    cell_data.method = 'TGFb_Junc';
    cell_data.jmag_max = 1e4;
end


% fname = ['sim2_cpm_tian_mcs',num2str(def.NRINC),...
%     '_',num2str(def.NVX),'x',num2str(def.NVY),...
%     '_scale',num2str(cellmod.params.model_params.scale),...
%     '.mat'];
%fname='test.mat';
%fname = ['test_' num2str(def.NVX) 'p_' num2str(def.NRINC) 'mcs_' num2str(tgfbval) 'tgfb_' num2str(tmaxval) 'tmax_' num2str(spreadval) 'spread_' 'div_' num2str(dividing) '_Jhalf' num2str(Jhalf) '_extraA' num2str(extrascale1) '_extraBT' num2str(extrascale2) '_ZEBs' num2str(pZEB) '_R200s' num2str(pR200) '_ZEB' num2str(intrazeb1) '_SNAIL' num2str(intrasnail2) '_kDAE' num2str(kdae_scale) '_' comment '.mat'];
%fname = ['test_' num2str(def.NVX) 'p_' num2str(def.NRINC) 'mcs_' num2str(tgfbval) 'tgfb_' num2str(spreadval) 'spread' num2str(dividing) 'div' '_extraA' num2str(extrascale1) '_extraBT' num2str(extrascale2) '_ZEB' num2str(intrazeb1) '_SNAIL' num2str(intrasnail2) '_kDAE.scale' num2str(kdae_scale) '_' comment '.mat'];
%fname = ['test_' num2str(def.NVX) 'p_' num2str(def.NRINC) 'mcs_' num2str(tgfbval) 'tgfb_' num2str(spreadval) 'spread' num2str(dividing) 'div' '_' comment '.mat'];
% dividing is now a vector
%fname = ['test_' num2str(def.NVX) 'p_' num2str(def.NRINC) 'mcs_' num2str(tgfbval) 'tgfb_' num2str(spreadval) 'spread_' comment '.mat'];
fname = ['test_' num2str(def.NVX) 'p_' num2str(def.NRINC) 'mcs_' num2str(tgfbval) 'tgfb_' comment '.mat'];


cpmfem_v2_1(fname, def, cellmod, phenoparams, cell_data);
pause(1);
%make_movie_or_snapshots(fname(1:end), 'movie',1:def.NRINC); %og was fname(1:end-4);
make_movie_or_snapshots_v3(fname(1:end), 'movie',1:def.NRINC); %og was fname(1:end-4);

toc

end