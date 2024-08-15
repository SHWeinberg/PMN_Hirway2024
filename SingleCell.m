%% TianXing BiosphyJ  2013 Model
% Shreyas Hirway
% 6/9/19

function [final] = SingleCell(init)
%% Initial Conditions
%x=[]; % T, s, S, R3, z, Z, R2, E, N Concentrations
tspan=init.time;% hours

%% Go through each cell

T00= init.T0;% T TGF-B
s0=init.s0;% s snail
S0=init.S0;% S SNAIL
R30=init.R30;% R3 miR-34
z0=init.z0;% z zeb
Z0=init.Z0;% Z ZEB
R20=init.R20;% R2 miR-200
E0=init.E0;% E E-cadherin
N0=init.N0;% N N-cadherin

x=zeros(1,9);
x(1)=T00;
x(2)=s0;
x(3)=S0;
x(4)=R30;
x(5)=z0;
x(6)=Z0;
x(7)=R20;
x(8)=E0;
x(9)=N0;
par=struct;
par.exo=init.exo;
par.Tbas=init.Tbas;
par.mir2=init.mir2;
% par.a=init.a;
% par.b=init.b;
% a = par.a; %lower bound for rand function
% b = par.b; %upper bound for rand function
par.T0=par.exo;%1; %Exogenous concentration of TGF-? 0?5 unit ******Can be changed to 3.0
par.k0T=par.Tbas;%k0T=0.06;  % Basal production rate of TGF-? 0.06 ?M/hr
par.kT=1.2;  % Production rate of TGF-? 1.2 ?M/hr
par.JT=0.06;  % Michaelis constant of miR200-dependent inhibition of TGF-? expression 0.06 ?M
par.kdT=0.6;  % Degradation rate of TGF-? 0.6/hr %Original in paper was 0.6/hr

par.k0s=0.0006;  % Basal transcription rate of snail1 0.0006 ?M/hr
par.ks=0.03;  %Transcription rate of snail1 0.03 ?M/hr
par.kds=0.09;  %Degradation rate of snail1 mRNA 0.09/hr
par.Js=1.6;  %Michaelis constant of TGF-?-dependent snail1 translation 1.6 ?M

par.kS=17;  %Translation rate of snail1 mRNA 17 ?M/hr
par.kdS=1.66;  % Degradation rate of SNAIL1 1.66/hr
par.JS=0.08;  % Michaelis constant of miR34-dependent inhibition of snail1 translation 0.08 ?M

par.k03=0.0012;  % Basal production rate of miR-34 0.0012 ?M/hr
par.k3=0.012;  % Production rate of miR-34 0.012 ?M/hr
par.kd3=0.035;  % Degradation rate of miR-34 0.035/hr
par.J13=0.15;  % Michaelis constant of SNAIL1-dependent inhibition of miR-34 production 0.15 ?M
par.J23=0.36;  % Michaelis constant of ZEB-dependent inhibition of miR-34 production 0.36 ?M

par.k0z=0.003;  % Basal transcription rate of zeb 0.003 ?M/hr
par.kz=0.06;  % Transcription rate of zeb 0.06 ?M/hr
par.Jz=3.5;  % Michaelis constant of SNAIL1-dependent zeb transcription 3.5 ?M
par.kdz=0.09;  % Degradation rate of zeb mRNA 0.09/hr

par.kdZ=1.66;  % Degradation rate of ZEB 1.66/hr
par.kZ=17;  % Translation rate of zeb mRNA 17 ?M/hr
par.JZ=0.06;  % Michaelis constant of miR34-dependent inhibition of zeb mRNA translation 0.06 ?M

par.k02=par.mir2;%k02=0.0002;  % Basal production rate of miR-200 0.0002 ?M/hr
par.k2=0.012;  % Production rate of miR-200 0.012 ?M/hr
par.kd2=0.035;  % Degradation rate of miR-200 0.035/hr Original: 0.035/hr
par.J12=5;  % Michaelis constant of SNAIL1-dependent inhibition of miR-200 production 5 ?M
par.J22=0.2;  % Michaelis constant of ZEB-dependent inhibition of miR-200 production 0.2 ?M

par.ke1=1;  % Production rate 1 of E-cadherin production 1 ?M/hr
par.ke2=0.6;  % Production rate 2 of E-cadherin production 0.6 ?M/hr
par.kde=0.5;  % Degradation rate of E-cadherin 0.5/hr
par.J1e=0.2;  % Michaelis constant of SNAIL1-dependent inhibition of E-cadherin production 0.2 ?M
par.J2e=0.5;  % Michaelis constant of ZEB-dependent inhibition of E-cadherin production 0.5 ?M

par.kn1=1;    % Production rate 1 of N-cadherin production 1 ?M/hr
par.kn2=0.6;  % Production rate 2 of N-cadherin production 0.6 ?M/hr
par.J1n=0.2;  % Michaelis constant of SNAIL1-dependent N-cadherin production 0.2?M
par.J2n=0.5;  % Michaelis constant of ZEB-dependent N-cadherin production 0.5 ?M
par.kdn=0.5;  % Degradation rate of N-cadherin 0.5/hr
par.scale=init.scale;
par.conv=init.conv;

%tgfb
% par.kT200;
% par.kbT;
% par.kE;
% par.kdE;
% par.kaNcad;
% par.kdAE;
% par.cell;
% Anon Function
%options=odeset('RelTol',1e-30,'AbsTol',1e-30);
[t,final]=ode15s(@(t,x) SingleCellFunc(tspan,x,par),tspan,x);

% Changed to 15s

end


