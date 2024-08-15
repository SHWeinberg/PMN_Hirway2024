kt2001=0.0943;%1.2;
kt2002=1.26;
kdt1=0.6;
kdt2=0.6;
kbt1=10.656;
kbt2=10.656;
kdae2=0.0576;
fdae=5;
kdae1=kdae2/fdae;
ka2=0.0833;
fa=10;
ka1=ka2/fa;
kde1=0.6;
kde2=0.6;
ke1=0.006;
ke2=0.0273;
exot1=0.16;
exot2=3;
sole1=0.01;
sole2=0.04;

a1=kdt1;
a2=kdt2;
b1=ka1*sole1-kt2001+kdae1*kdt1/kbt1+exot1*kdt1;
b2=ka2*sole2-kt2002+kdae2*kdt2/kbt2+exot2*kdt2;
c1=-kdae1*kt2001/kbt1-exot1*kt2001;
c2=-kdae2*kt2002/kbt2-exot2*kt2002;
d1=b1^2-4*a1*c1;
d2=b2^2-4*a2*c2;

nump1=-b1+sqrt(d1);
nump2=-b2+sqrt(d2);
numm1=-b1-sqrt(d1);
numm2=-b2-sqrt(d2);
den1=2*a1;
den2=2*a2;
rootp1=nump1/den1;
rootp2=nump2/den2;
rootm1=numm1/den1;
rootm2=numm2/den2;
solt1=rootp1;
solt2=rootp2;
asme1=(ka1*sole1)/(kdae1+kbt1*(solt1+exot1));
asme2=(ka2*sole2)/(kdae2+kbt2*(solt2+exot2));

ecmt1=(kbt1*(solt1+exot1)*asme1)/kdae1;
ecmt2=(kbt2*(solt2+exot2)*asme2)/kdae2;

%% Using functions Low Ncad
syms solT asmE ecmT
eq1= solT - (kt2001/(kdt1+kbt1*asmE));
eq2= asmE-(ka1*sole1)/(kdae1+kbt1*(solT+exot1));
eq3=ecmT-(kbt1*(solT+exot1)*asmE)/kdae1;
sol = solve(eq1,eq2,eq3);
disp('SolT is ')
double(sol.solT)
disp('asmE is ')
double(sol.asmE)
disp('ecmT is ')
double(sol.ecmT)

%% Using functions High Ncad
syms solT asmE ecmT
eq1= solT - (kt2002/(kdt2+kbt2*asmE));
eq2= asmE-(ka2*sole2)/(kdae2+kbt2*(solT+exot2));
eq3=ecmT-(kbt2*(solT+exot2)*asmE)/kdae2;
sol = solve(eq1,eq2,eq3);
disp('SolT is ')
double(sol.solT)
disp('asmE is ')
double(sol.asmE)
disp('ecmT is ')
double(sol.ecmT)