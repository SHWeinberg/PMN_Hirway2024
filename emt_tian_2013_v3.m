function dX = emt_tian_2013_v3(t,X,param,cellparam)
% to match update_cellmarkers_v2 file
% params - structure of parameters that does NOT vary b/w cells
% cellparam - structure of parameters that do vary b/w cells

% Units:
% concentration - uM
% time - min (parameters are in hr, ODE converted into min at end of
% function)
%
% param - structure of model parameters, same for all cells, constant in
% time
% cellparam - structure of model parameters that may differ between cells
% and/or vary in time (ex. exog TGFb)

T = X(1);
s = X(2);
S = X(3);
R3 = X(4);
z = X(5);
Z = X(6);
R2 = X(7);
E = X(8);
N = X(9);

TGF0 = cellparam.TGF0;
J=cellparam.J;

% parameters
k0T = param.k0T;
kT = param.kT;
JT = param.JT;   
kdT = param.kdT; 

k0s = param.k0s;
ks = param.ks;
Js = param.Js;
kds = param.kds;

kS = param.kS;
JS = param.JS;
kdS = param.kdS;

k03 = param.k03;
k3 = param.k3;
J13 = param.J13;
J23 = param.J23;
kd3 = param.kd3;

k0z = param.k0z;
kz = param.kz;
Jz = param.Jz;
kdz = param.kdz;

kZ = param.kZ;
JZ = param.JZ;
kdZ = param.kdZ;

k02 = param.k02;
k2 = param.k2;
J12 = param.J12;
J22 = param.J22;
kd2 = param.kd2;

ke1 = param.ke1;
J1e = param.J1e;
ke2 = param.ke2;
J2e = param.J2e;
kde = param.kde;

kn1 = param.kn1;
J1n = param.J1n;
kn2 = param.kn2;
J2n = param.J2n;
kdn = param.kdn;

nr2 = param.nr2;
nt = param.nt;
ns = param.ns;
nz = param.nz;
nr3 = param.nr3;

Jhalf=param.Jhalf;
nf=param.nf;
nf2=param.nf2;
fmax=param.fmax;

f=fmax/(1+(J/Jhalf)^nf2);
dT = k0T + kT/(1+(R2/JT)^nr2) - kdT*T;
ds = k0s + ks*(((T+TGF0)/Js)^nt/(1+((T+TGF0)/Js)^nt) + (f^nf)/(1+(f^nf))) - kds*s; 
dS = kS*s/(1+(R3/JS)^nr3) - kdS*S;
dR3 = k03 + k3/(1+(S/J13)^ns + (Z/J23)^nz) - kd3*R3;
dz = k0z + kz*(S/Jz)^ns/(1+(S/Jz)^ns) - kdz*z;
dZ = kZ*z/(1+(R2/JZ)^nr2) - kdZ*Z;
dR2 = k02 + k2/(1+(S/J12)^ns + (Z/J22)^nz) - kd2*R2;
dE = ke1/(1+(S/J1e)^ns) + ke2/(1+(Z/J2e)^nz) - kde*E;
dN = kn1*(S/J1n)^ns/(1+(S/J1n)^ns) + kn2*(Z/J2n)^nz/(1+(Z/J2n)^nz) - kdn*N;


dX = param.scale*[dT; ds; dS; dR3; dz; dZ; dR2; dE; dN]/60; % uM/min