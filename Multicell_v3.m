%%TianXing BiosphyJ  2013 Model
%%Shreyas Hirway
%Modified 10/26

function [final,extra] = Multicell_v3(init,cells)
%% Initial Conditions
%x=[]; % T, s, S, R3, z, Z, R2, E, N Concentrations
tspan=init.time;% hours

%% Go through each cell
par=struct;
par.mask=init.mask;
par.exo=init.exo;
par.scale=init.scale;
par.conv=init.conv;
%par.Tbas=init.Tbas;
%par.mir2=init.mir2;
par.cells=cells;
par.nx=size(par.mask,1)^2;
nx=par.nx;
%par.n=par.nx^2+cells;

%T00=zeros(cells,1);
s0=zeros(cells,1);
S0=zeros(cells,1);
R30=zeros(cells,1);
z0=zeros(cells,1);
Z0=zeros(cells,1);
R20=zeros(cells,1);
E0=zeros(cells,1);
N0=zeros(cells,1);


solT=zeros(nx,1);
solE=zeros(nx,1);
asmE=zeros(nx,1);
ecmT=zeros(nx,1);
totalE=zeros(nx,1);
totalT=zeros(nx,1);

%T00(:)= init.T0;% T TGF-B
s0(:)=init.s0;% s snail
S0(:)=init.S0;% S SNAIL
R30(:)=init.R30;% R3 miR-34
z0(:)=init.z0;% z zeb
Z0(:)=init.Z0;% Z ZEB
R20(:)=init.R20;% R2 miR-200
E0(:)=init.E0;% E E-cadherin
N0(:)=init.N0;% N N-cadherin

solT(:)=init.solT; % soluble TgfB
solE(:)=init.solE; % soluble ECM
asmE(:)=init.asmE; % assembled ECM
ecmT(:)=init.ecmT; % ECM bound TgfB

totalE(:)=asmE(:)+ecmT(:); % new Total ECM 
totalT(:)=solT(:)+ecmT(:); % new Total TGFB

par.param_js=init.param_js;
par.param_ke=init.param_ke;
par.param_ka=init.param_ka;
par.param_kdae=init.param_kdae;
% par.partialscaleZEB=init.partialscaleZEB;
% par.partialscalemir200=init.partialscalemir200;

par.parametermatrix=init.parameters;
par.extra=init.extra;
% par.cell;

% Creating Input Vector
x=[];
%x=[ solT; solE; asmE; ecmT; s0; S0; R30; z0; Z0; R20; E0; N0];

% Input vector with new extracellular terms
x=[solE; totalE; totalT; s0; S0; R30; z0; Z0; R20; E0; N0];

% Anon Function

[t,final]=ode15s(@(t,x) TestCellFunc_v4(tspan,x,par),tspan,x);

% Calculate solT, asmE and ecmT here
totalEv=final(:,1*nx+1:2*nx);
totalTv=final(:,2*nx+1:3*nx);
Ncad_val=zeros(size(final,1),nx);
NVal=final(:,3*nx+par.cells*7+1:end);
for i=1:par.cells
    Ncad_temp=NVal(:,i);
    for time=1:size(final,1)
    Ncad_val(time,find(par.mask==i))=Ncad_temp(time);
    end
end
%Ncad_val=reshape(Ncad_val,nx,1);
kbT=par.extra.kbT;
%kdAE=par.extra.kdAE;
kdAE2=par.extra.kdAE2;%(720/500/25)*param_kdae;
kdAE1=par.extra.kdAE1;%kdAE2/fdAE; % from ECM degradation graph
kdAE=(kdAE2-kdAE1).*Ncad_val/3.1515 + kdAE1; 
extra=[];
ecmTv=(kbT.*(totalTv+totalEv)+kdAE - sqrt( (kbT.*(totalTv+totalEv)+kdAE).^2 - 4.*kbT.^2.*totalTv.*totalEv))./(2*kbT);
solTv=totalTv-ecmTv;
asmEv=totalEv-ecmTv;
extra=[solTv,asmEv,ecmTv];
% Changed to 15s
end