function ddt = SingleCellFunc(t,mat,par)


x=mat;

T0=par.T0;%1; %Exogenous concentration of TGF-? 0?5 unit ******Can be changed to 3.0
k0T=par.k0T;%k0T=0.06;  % Basal production rate of TGF-? 0.06 ?M/hr
kT=par.kT;  % Production rate of TGF-? 1.2 ?M/hr
JT=par.JT;  % Michaelis constant of miR200-dependent inhibition of TGF-? expression 0.06 ?M
kdT=par.kdT;  % Degradation rate of TGF-? 0.6/hr %Original in paper was 0.6/hr

k0s=par.k0s;  % Basal transcription rate of snail1 0.0006 ?M/hr
ks=par.ks;  %Transcription rate of snail1 0.03 ?M/hr
kds=par.kds;  %Degradation rate of snail1 mRNA 0.09/hr
Js=par.Js;  %Michaelis constant of TGF-?-dependent snail1 translation 1.6 ?M

kS=par.kS;  %Translation rate of snail1 mRNA 17 ?M/hr
kdS=par.kdS;  % Degradation rate of SNAIL1 1.66/hr
JS=par.JS;  % Michaelis constant of miR34-dependent inhibition of snail1 translation 0.08 ?M

k03=par.k03;  % Basal production rate of miR-34 0.0012 ?M/hr
k3=par.k3;  % Production rate of miR-34 0.012 ?M/hr
kd3=par.kd3;  % Degradation rate of miR-34 0.035/hr
J13=par.J13;  % Michaelis constant of SNAIL1-dependent inhibition of miR-34 production 0.15 ?M
J23=par.J23;  % Michaelis constant of ZEB-dependent inhibition of miR-34 production 0.36 ?M

k0z=par.k0z;  % Basal transcription rate of zeb 0.003 ?M/hr
kz=par.kz;  % Transcription rate of zeb 0.06 ?M/hr
Jz=par.Jz;  % Michaelis constant of SNAIL1-dependent zeb transcription 3.5 ?M
kdz=par.kdz;  % Degradation rate of zeb mRNA 0.09/hr

kdZ=par.kdZ;  % Degradation rate of ZEB 1.66/hr
kZ=par.kZ;  % Translation rate of zeb mRNA 17 ?M/hr
JZ=par.JZ;  % Michaelis constant of miR34-dependent inhibition of zeb mRNA translation 0.06 ?M

k02=par.k02;%k02=0.0002;  % Basal production rate of miR-200 0.0002 ?M/hr
k2=par.k2;  % Production rate of miR-200 0.012 ?M/hr
kd2=par.kd2;  % Degradation rate of miR-200 0.035/hr Original: 0.035/hr
J12=par.J12;  % Michaelis constant of SNAIL1-dependent inhibition of miR-200 production 5 ?M
J22=par.J22;  % Michaelis constant of ZEB-dependent inhibition of miR-200 production 0.2 ?M

ke1=par.ke1;  % Production rate 1 of E-cadherin production 1 ?M/hr
ke2=par.ke2;  % Production rate 2 of E-cadherin production 0.6 ?M/hr
kde=par.kde;  % Degradation rate of E-cadherin 0.5/hr
J1e=par.J1e;  % Michaelis constant of SNAIL1-dependent inhibition of E-cadherin production 0.2 ?M
J2e=par.J2e;  % Michaelis constant of ZEB-dependent inhibition of E-cadherin production 0.5 ?M

kn1=par.kn1;    % Production rate 1 of N-cadherin production 1 ?M/hr
kn2=par.kn2;  % Production rate 2 of N-cadherin production 0.6 ?M/hr
J1n=par.J1n;  % Michaelis constant of SNAIL1-dependent N-cadherin production 0.2?M
J2n=par.J2n;  % Michaelis constant of ZEB-dependent N-cadherin production 0.5 ?M
kdn=par.kdn;  % Degradation rate of N-cadherin 0.5/hr

nt=2;  % Hill coefficient of TGF-?-dependent SNAIL1 expression 2
nr2=2;  % Hill coefficient of miR-200-dependent inhibition 2
nr3=2;  % Hill coefficient of miR-34-dependent inhibition 2
nz=2;  % Hill coefficient of ZEB-dependent inhibition 2
ns=2;  % Hill coefficient of SNAIL1-dependent activation or inhibition 2

%TGFB
% kT200=par.kT200;
% kbT=par.kbT;
% kE=par.kE;
% kdE=par.kdE;
% kaNcad=par.kaNcad;
% kdAE=par.kdAE;
% cell=par.cell;
%x=[]; % T, s, S, R3, z, Z, R2, E, N Concentrations

%% Equations


A= (kT ./ (1+ (x(7)/JT).^nr2));

Eq1 = (k0T+A-kdT*x(1));
B= (((x(1)+ T0)./Js).^nt) ./ (1+(((x(1)+T0)./Js).^nt));
Eq2=  (k0s+ks*B-kds*x(2));
C=1./ (1+ (x(4)./JS).^nr3);
Eq3=  (kS*x(2).*C-kdS*x(3));
D= 1./ ( 1+ ((x(3)./J13).^ns) + ((x(6)./J23).^nz));
Eq4= (k03+k3*D-kd3*x(4));
E=((x(3)./Jz).^ns) ./ (1+((x(3)./Jz).^ns));
Eq5= (k0z+kz*E-kdz*x(5));
F=1./(1+((x(7)./JZ).^nr2));
Eq6= (kZ*x(5).*F-kdZ*x(6));
G=1./ ( 1+ ((x(3)./J12).^ns) + ((x(6)./J22).^nz));
Eq7= (k02+k2*G-kd2*x(7));
H=1./(1+((x(3)./J1e).^ns));
I=1./(1+((x(6)./J2e).^nz));
Eq8= (ke1*H+ke2*I-kde*x(8));
J=((x(3)./J1n).^ns) ./ (1+((x(3)./J1n).^ns));
K=((x(6)./J2n).^nz) ./ (1+((x(6)./J2n).^nz));
Eq9= (kn1*J+kn2*K-kdn*x(9));

% Tgfb dynamics
% Eq10=kT200*(1-kroneckerDelta(cell))-kdT*x(10)-kbT*x(10)*x(12);
% Eq11=kE*(1-kroneckerDelta(cell))-kdE*x(11)-kaNcad*(1-kroneckerDelta(cell))*x(11);
% Eq12=kaNcad*(1-kroneckerDelta(cell))*x(11)-kdAE*x(12)-kbT*x(10)*x(12);
% Eq13=kbT*x(10)*x(12)-kdE*x(13);

scaleF=par.scale; %Scaling factor to make timeline similar to in vitro
conv=par.conv; %Conversion factor for time
ddt=conv*scaleF*[Eq1 Eq2 Eq3 Eq4 Eq5 Eq6 Eq7 Eq8 Eq9]';


end