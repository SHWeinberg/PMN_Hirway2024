function ddt = CellFunc_v5(t,mat,par)
x=mat;%reshape(mat,par.n,13);
param_js=par.extra.param_js;
param_ke=par.extra.param_ke;
param_ka=par.extra.param_ka;
param_kdae=par.extra.param_kdae;
% partialscaleZEB=par.partialscaleZEB;
% partialscalemir200=par.partialscalemir200;

nx=par.nx;
cells=par.cells;
%mask=zeros(size(par.mask,1)^2,1);
%maskC=[mask;ones(par.cells,1)];
%maskG=reshape(par.mask,size(par.mask,1)^2,1); inputted mask will already
%be a vector
maskG=par.mask;

T0=zeros(sqrt(nx),sqrt(nx));
%T0(35:65,35:65)=par.cellTgfb(1); %ONLY CENTER of grid gets tgfb
T0(:)=par.cellTgfb(1);%1; %Exogenous concentration of TGF-? 0?5 unit ******Can be changed to 3.0
% Since every cell gets the same value, T0(:)=just 1 value. can later
% change to a matrix
T0=reshape(T0,nx,1);

%T0 is exogenous TGFB

% Intracellular parameters- new order compatible with Mario's parameters
% and tian_baseline_v3

k0T=par.parametermatrix(1,:)';%0.06;  % Basal production rate of TGF-? 0.06 ?M/hr
k0s=par.parametermatrix(2,:)';%0.0006;  % Basal transcription rate of snail1 0.0006 ?M/hr
k03=par.parametermatrix(3,:)';%0.0012;  % Basal production rate of miR-34 0.0012 ?M/hr
k0z=par.parametermatrix(4,:)';%0.003;  % Basal transcription rate of zeb 0.003 ?M/hr
k02=par.parametermatrix(5,:)';%0.0002;  % Basal production rate of miR-200 0.0002 ?M/hr

kT=par.parametermatrix(6,:)';%1.2;  % Production rate of TGF-? 1.2 ?M/hr
ks=par.parametermatrix(7,:)';%0.03;  %Transcription rate of snail1 0.03 ?M/hr
kS=par.parametermatrix(8,:)';%17;  %Translation rate of snail1 mRNA 17 ?M/hr
k3=par.parametermatrix(9,:)';%0.012;  % Production rate of miR-34 0.012 ?M/hr
kz=par.parametermatrix(10,:)';%0.06;  % Transcription rate of zeb 0.06 ?M/hr
kZ=par.parametermatrix(11,:)';%17*partialscaleZEB;  % Translation rate of zeb mRNA 17 ?M/hr *Could decrease this for partial state- SUH 111020
k2=par.parametermatrix(12,:)';%0.012*partialscalemir200;  % Production rate of miR-200 0.012 ?M/hr *Could increase this for partial state- SUH 111020

kdT=par.parametermatrix(13,:)';%0.6;  % Degradation rate of TGF-? 0.6/hr
kds=par.parametermatrix(14,:)';%0.09;*intrasnail2;   %Degradation rate of snail1 mRNA 0.09/hr
kdS=par.parametermatrix(15,:)';%1.66;  % Degradation rate of SNAIL1 1.66/hr
kdz=par.parametermatrix(16,:)';%0.09;*intrazeb1;   % Degradation rate of zeb mRNA 0.09/hr
kdZ=par.parametermatrix(17,:)';%1.66;  % Degradation rate of ZEB 1.66/hr
kd3=par.parametermatrix(18,:)';%0.035;  % Degradation rate of miR-34 0.035/hr
kd2=par.parametermatrix(19,:)';%0.035;  % Degradation rate of miR-200 0.035/hr

JT=par.parametermatrix(20,:)';%0.06;  % Michaelis constant of miR200-dependent inhibition of TGF-? expression 0.06 ?M
Js=par.parametermatrix(21,:)';%1.6;  %Michaelis constant of TGF-?-dependent snail1 translation 1.6 ?M
JS=par.parametermatrix(22,:)';%0.08;  % Michaelis constant of miR34-dependent inhibition of snail1 translation 0.08 ?M
J13=par.parametermatrix(23,:)';%0.15;  % Michaelis constant of SNAIL1-dependent inhibition of miR-34 production 0.15 ?M
J23=par.parametermatrix(24,:)';%0.36;  % Michaelis constant of ZEB-dependent inhibition of miR-34 production 0.36 ?M
Jz=par.parametermatrix(25,:)';%3.5;  % Michaelis constant of SNAIL1-dependent zeb transcription 3.5 ?M
JZ=par.parametermatrix(26,:)';%0.06;  % Michaelis constant of miR34-dependent inhibition of zeb mRNA translation 0.06 ?M
J12=par.parametermatrix(27,:)';%5;  % Michaelis constant of SNAIL1-dependent inhibition of miR-200 production 5 ?M
J22=par.parametermatrix(28,:)';%0.2;  % Michaelis constant of ZEB-dependent inhibition of miR-200 production 0.2 ?M

ke1=par.parametermatrix(29,:)';%1;  % Production rate 1 of E-cadherin production 1 ?M/hr
J1e=par.parametermatrix(30,:)';%0.2;  % Michaelis constant of SNAIL1-dependent inhibition of E-cadherin production 0.2 ?M
ke2=par.parametermatrix(31,:)';%0.6;  % Production rate 2 of E-cadherin production 0.6 ?M/hr
J2e=par.parametermatrix(32,:)';%0.5;  % Michaelis constant of ZEB-dependent inhibition of E-cadherin production 0.5 ?M
kde=par.parametermatrix(33,:)';%0.5;  % Degradation rate of E-cadherin 0.5/hr

kn1=par.parametermatrix(34,:)';%1;  % Production rate 1 of N-cadherin production 1 ?M/hr
J1n=par.parametermatrix(35,:)';%0.2;  % Michaelis constant of SNAIL1-dependent N-cadherin production 0.2?M
kn2=par.parametermatrix(36,:)';%0.6;  % Production rate 2 of N-cadherin production 0.6 ?M/hr
J2n=par.parametermatrix(37,:)';%0.5;  % Michaelis constant of ZEB-dependent N-cadherin production 0.5 ?M
kdn=par.parametermatrix(38,:)';%0.5;  % Degradation rate of N-cadherin 0.5/hr

nt=par.parametermatrix(39,:)';%2;  % Hill coefficient of TGF-?-dependent SNAIL1 expression 2
ns=par.parametermatrix(40,:)';%2;  % Hill coefficient of SNAIL1-dependent activation or inhibition 2
nz=par.parametermatrix(41,:)';%2;  % Hill coefficient of ZEB-dependent inhibition 2
nr2=par.parametermatrix(42,:)';%2;  % Hill coefficient of miR-200-dependent inhibition 2
nr3=par.parametermatrix(43,:)';%2;  % Hill coefficient of miR-34-dependent inhibition 2




% par.kbT=1;
% par.kE=0.02; % from 20nm ECM fibronectin prod
% par.kdE=1.2;
% par.kaNcad=1;
% par.kdAE=1.2;


%Junctional Forces related stuff
cellJ=par.cellJ;
Jhalf=par.extra.Jhalf;
nf=par.extra.nf;
nf2=par.extra.nf2;
fmax=par.extra.tmax./Js;

%% Equations
%disp('Running...')
% solTv=x(1:nx);
% solEv=x(nx+1:2*nx);
% asmEv=x(2*nx+1:3*nx);
% ecmTv=x(3*nx+1:4*nx);
% sVal=x(4*nx+1:4*nx+1*cells);
% SVal=x(4*nx+1*cells+1: 4*nx+2*cells);
% R3Val=x(4*nx+2*cells+1: 4*nx+3*cells);
% zVal=x(4*nx+3*cells+1: 4*nx+4*cells);
% ZVal=x(4*nx+4*cells+1: 4*nx+5*cells);
% R2Val=x(4*nx+5*cells+1: 4*nx+6*cells);
% EVal=x(4*nx+6*cells+1: 4*nx+7*cells);
% NVal=x(4*nx+7*cells+1: 4*nx+8*cells);

solEv=x(1:nx);
totalEv=x(1*nx+1:2*nx);
totalTv=x(2*nx+1:3*nx);

sVal=x(3*nx+1:3*nx+1*cells);
SVal=x(3*nx+1*cells+1: 3*nx+2*cells);
R3Val=x(3*nx+2*cells+1: 3*nx+3*cells);
zVal=x(3*nx+3*cells+1: 3*nx+4*cells);
ZVal=x(3*nx+4*cells+1: 3*nx+5*cells);
R2Val=x(3*nx+5*cells+1: 3*nx+6*cells);
EVal=x(3*nx+6*cells+1: 3*nx+7*cells);
NVal=x(3*nx+7*cells+1: 3*nx+8*cells);


Ncad_val=zeros(sqrt(nx),sqrt(nx));
for i=1:cells
    Ncad_temp=NVal(i);
    Ncad_val(find(maskG==i))=Ncad_temp;
end
Ncad_val=reshape(Ncad_val,nx,1);

% Previous Version with scaling extracellular rates- SUH 10/25/21
extrascaling1=par.extra.scaling1;
extrascaling2=par.extra.scaling2;
kdae_scale=par.extra.kdae_scale;
ke_scale=par.extra.ke_scale;
%extracellular dynamics TGFB
kbT=par.extra.kbT*extrascaling2;%*extrascaling;%10.656; % 0.008- calculated from 12th-14th paper
kE1=par.extra.kE1;%0.006; % from 10-40nm ECM fibronectin prod
kE2=par.extra.kE2;%0.0273;
kE=((kE2-kE1).*Ncad_val/3.1515 + kE1)*param_ke; % Uses Ncad so calculated here
kE=kE*ke_scale;
kdE=par.extra.kdE;%0.6; %from Tian Model
fa=par.extra.fa;%10; % Scaling Factor
ka1=par.extra.ka1;%1/(fa*12);
ka2=par.extra.ka2;%(1/12)*param_ka;
ka=(ka2-ka1).*Ncad_val/3.1515 + ka1; % Uses Ncad so calculated here
ka=ka*extrascaling1;
fdAE=par.extra.fdAE;%5; % Scaling Factor
kdAE2=par.extra.kdAE2;%(720/500/25)*param_kdae;
kdAE1=par.extra.kdAE1;%kdAE2/fdAE; % from ECM degradation graph
kdAE=(kdAE2-kdAE1).*Ncad_val/3.1515 + kdAE1; % Uses Ncad so calculated here
kdAE=kdAE*kdae_scale;
% ****************************


% ********************MAJOR ISSUE WITH vectorized extracellular rates
% causing cells to freeze- SUH 10/25/21
extrascaling1=par.extra.scaling1;
extrascaling2=par.extra.scaling2;
kdae_scale=par.extra.kdae_scale;
ke_scale=par.extra.ke_scale;
%extracellular dynamics TGFB
kbT=par.parametermatrix(44,:)'*extrascaling2;%*extrascaling;%10.656; % 0.008- calculated from 12th-14th paper
kE1=par.parametermatrix(45,:)';%0.006; % from 10-40nm ECM fibronectin prod
kE2=par.parametermatrix(46,:)';%0.0273;
kdE=par.parametermatrix(47,:)';%0.6; %from Tian Model
fa=par.parametermatrix(48,:)';%10; % Scaling Factor
ka1=par.parametermatrix(49,:)';%1/(fa*12);
ka2=par.parametermatrix(50,:)';%(1/12)*param_ka;
fdAE=par.parametermatrix(51,:)';%5; % Scaling Factor
kdAE2=par.parametermatrix(52,:)';%(720/500/25)*param_kdae;
kdAE1=par.parametermatrix(53,:)';%kdAE2/fdAE; % from ECM degradation graph

kbT_mask=zeros(sqrt(nx),sqrt(nx));
kbT_mask(:)=par.extra.kbT*extrascaling2;
kdE_mask=zeros(sqrt(nx),sqrt(nx));
kdE_mask(:)=par.extra.kdE;
kE1_mask=zeros(sqrt(nx),sqrt(nx));
kE1_mask(:)=par.extra.kE1;
kE2_mask=zeros(sqrt(nx),sqrt(nx));
ka1_mask=zeros(sqrt(nx),sqrt(nx));
ka1_mask(:)=par.extra.ka1;
ka2_mask=zeros(sqrt(nx),sqrt(nx));
kdAE1_mask=zeros(sqrt(nx),sqrt(nx));
kdAE1_mask(:)=par.extra.kdAE1;
kdAE2_mask=zeros(sqrt(nx),sqrt(nx));
for i=1:cells
    kbT_temp=kbT(i);
    kbT_mask(find(maskG==i))=kbT_temp;
    
    kdE_temp=kdE(i);
    kdE_mask(find(maskG==i))=kdE_temp;
    
    kE1_temp=kE1(i);
    kE1_mask(find(maskG==i))=kE1_temp;
    
    kE2_temp=kE2(i);
    kE2_mask(find(maskG==i))=kE2_temp;
    
    ka1_temp=ka1(i);
    ka1_mask(find(maskG==i))=ka1_temp;
    
    ka2_temp=ka2(i);
    ka2_mask(find(maskG==i))=ka2_temp;
    
    kdAE1_temp=kdAE1(i);
    kdAE1_mask(find(maskG==i))=kdAE1_temp;
    
    kdAE2_temp=kdAE2(i);
    kdAE2_mask(find(maskG==i))=kdAE2_temp;
end

kbT_mask=reshape(kbT_mask,nx,1);
kdE_mask=reshape(kdE_mask,nx,1);
kE1_mask=reshape(kE1_mask,nx,1);
kE2_mask=reshape(kE2_mask,nx,1);
ka1_mask=reshape(ka1_mask,nx,1);
ka2_mask=reshape(ka2_mask,nx,1);
kdAE1_mask=reshape(kdAE1_mask,nx,1);
kdAE2_mask=reshape(kdAE2_mask,nx,1);


kbT=kbT_mask;
kdE=kdE_mask;
kE=((kE2_mask-kE1_mask).*Ncad_val/3.1515 + kE1_mask)*param_ke; % Uses Ncad so calculated here
kE=kE*ke_scale;

ka=(ka2_mask-ka1_mask).*Ncad_val/3.1515 + ka1_mask; % Uses Ncad so calculated here
ka=ka*extrascaling1;

kdAE=(kdAE2_mask-kdAE1_mask).*Ncad_val/3.1515 + kdAE1_mask; % Uses Ncad so calculated here
kdAE=kdAE*kdae_scale;
%disp([num2str(kdae_scale)]);
%kdAE=kdAE*1.1; % Need to figure out degradation rate ******************************* SUH 033121

%********************************* New version with vectorized
%extracellular rates ends here

%kdAE_dup=(kdAE2-kdAE1) + kdAE1;
MR200_val=zeros(sqrt(nx),sqrt(nx));
for i=1:cells
    M200_temp=R2Val(i);
    MR200_val(find(maskG==i))=M200_temp;
end
MR200_val=reshape(MR200_val,nx,1);

%grid matrices for k0T, kT, JT, nr2
k0T_val=zeros(sqrt(nx),sqrt(nx));
kT_val=zeros(sqrt(nx),sqrt(nx));
JT_val=zeros(sqrt(nx),sqrt(nx));
nr2_val=zeros(sqrt(nx),sqrt(nx));
kdT_val=zeros(sqrt(nx),sqrt(nx));
for i=1:cells
    k0T_temp=k0T(i);
    k0T_val(find(maskG==i))=k0T_temp;

    kT_temp=kT(i);
    kT_val(find(maskG==i))=kT_temp;
    
    JT_temp=JT(i);
    JT_val(find(maskG==i))=JT_temp; 
    
    nr2_temp=nr2(i); %Was accidentally k0T
    nr2_val(find(maskG==i))=nr2_temp;  
    
    
end

kdT_temp=kdT(i);
    kdT_val(:)=kdT_temp; % Changed to not be cell dependent- uniform degradation everywhere-SUH 111220

k0T_val=reshape(k0T_val,nx,1);
kT_val=reshape(kT_val,nx,1);
JT_val=reshape(JT_val,nx,1);
nr2_val=reshape(nr2_val,nx,1);

kdT_val=reshape(kdT_val,nx,1);



kT200 = k0T_val + kT_val./(1 + ([MR200_val]./JT_val).^nr2_val);
kT200(isnan(kT200))=0;

%Extracellular
%dsolT=kT200.*(sign(maskG))-kdT_val.*solTv-kbT.*solTv.*asmEv;
dsolE=kE.*(sign(maskG))-kdE.*solEv-ka.*(sign(maskG)).*solEv;
%dasmeE=ka.*(sign(maskG)).*solEv-kdAE.*(sign(maskG)).*asmEv-kbT.*(solTv+T0).*asmEv;
%decmT=kbT.*(solTv+T0).*asmEv-kdAE.*ecmTv;
dtotalE=ka.*(sign(maskG)).*solEv-kdAE.*(sign(maskG)).*totalEv; % sign(maskG) added to kdAE SUH 7/16/21
% Quadratic Formula here ***************************************
ecmTv=(kbT.*(totalTv+totalEv)+kdAE.*(sign(maskG)) - sqrt( (kbT.*(totalTv+totalEv)+kdAE.*(sign(maskG))).^2 - 4.*kbT.^2.*totalTv.*totalEv))./(2*kbT); % added .*(sign(maskG)) to kdAE- SUH 10/25/21
solTv=totalTv-ecmTv;
dtotalT=kT200 + T0.*kdT_val - kdT_val.*(sign(maskG)).*solTv-kdAE.*(sign(maskG)).*ecmTv; % kt,exo changed from T0/kdt to T0*kdt
% sign(maskG) added to kdt_val and kdAE SUH 7/16/21

% dsolT=dsolT*extrascaling; 
% dsolE=dsolE*extrascaling;
% dasmeE=dasmeE*extrascaling;
% decmT=decmT*extrascaling;
% Ncad_val=zeros(sqrt(nx),sqrt(nx));
% for i=1:cells
%     Ncad_temp=NVal(i);
%     Ncad_val(find(maskG==i))=Ncad_temp;
% end
% Ncad_val=reshape(Ncad_val,nx,1);

Ecmval=[];
for i=1:cells
    Gval=find(maskG==i);
    Ecmval=[Ecmval; mean(ecmTv(Gval))];
end

% Ecmval(:)=0; %Testing without EcmT
%Intracellular
f=fmax./(1+(cellJ/Jhalf).^nf2);
%f(:)=0; % Testing without Jforces

Jss=(Js/(86.4865))*param_js; %Modified Js value- 86.4865 added as scaling factor SUH:5/28/20
B= (((Ecmval)./Jss).^nt) ./ (1+(((Ecmval)./Jss).^nt)) + (f.^nf)./(1+(f.^nf)); %86.4865 added as scaling factor SUH:5/21/20
ds= (k0s+ks.*B-kds.*sVal);
C=1./ (1+ (R3Val./JS).^nr3);
dS= (kS.*sVal.*C-kdS.*SVal);
D= 1./ ( 1+ ((SVal./J13).^ns) + ((ZVal./J23).^nz));
dR3=(k03+k3.*D-kd3.*R3Val);
E=((SVal./Jz).^ns) ./ (1+((SVal./Jz).^ns));
dz=(k0z+kz.*E-kdz.*zVal);
F=1./(1+((R2Val./JZ).^nr2));
dZ=(kZ.*zVal.*F-kdZ.*ZVal);
G=1./ ( 1+ ((SVal./J12).^ns) + ((ZVal./J22).^nz));
dR2=(k02+k2.*G-kd2.*R2Val);
H=1./(1+((SVal./J1e).^ns));
I=1./(1+((ZVal./J2e).^nz));
dE=(ke1.*H+ke2.*I-kde.*EVal);
J=((SVal./J1n).^ns) ./ (1+((SVal./J1n).^ns));
K=((ZVal./J2n).^nz) ./ (1+((ZVal./J2n).^nz));
dN=(kn1.*J+kn2.*K-kdn.*NVal);
%disp(['Time:' num2str(t)]);

scale=par.scale; %Scaling factor to make timeline similar to in vitro
conv=par.conv; %Conversion factor for time
% ddt=conv*scale*[dsolT; dsolE; dasmeE; decmT; ds; dS; dR3; dz; dZ; dR2; dE; dN];

% New Equation with Total ECM & TGFB
ddt=conv*scale*[dsolE; dtotalE; dtotalT; ds; dS; dR3; dz; dZ; dR2; dE; dN];
%ddt=ddt';

end