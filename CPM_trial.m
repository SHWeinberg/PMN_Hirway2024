dimval=100;
mcsval=252;
tgfbval=0;
spreadval=0.015*4;
tmaxval=0;
Jmin=4819*13;
dividing=1.25;%/4;
al=0.05;
n=2;
Jhalf0=Jmin/((al^-1-1)^(1/n));
Tmax=2*.15;
tmax0=1*.15;
Jmin= Jhalf0*(al^-1 -1) ^(1/n) / (2*Tmax/tmax0 - 1)^(1/n);
Jhalf=Jmin/((al^-1-1)^(1/n));
comment='';
pZEB=0.8;
pR200=1.2;
tmaxval=Tmax;
extrascale=1;
%%
clc
run_single_cpm_fem_ode(dimval,mcsval,tgfbval,spreadval,tmaxval,Jhalf,dividing,extrascale,pZEB,pR200,comment);
pause(1)
close all

%%

for sc=[0.88]
    for tgfbval=[0 50]
clc
extrascale=sc;
run_single_cpm_fem_ode(dimval,mcsval,tgfbval,spreadval,tmaxval,Jhalf,dividing,extrascale,pZEB,pR200,comment);
pause(1)
close all
    end
end

%%
tic
for i=[0, 1.5, 3]
    dividing=1/4;
dimval=100;
mcsval=2000;
spreadval=0.15;
tmaxval=1.5;
Jhalf= 638.2921;
run_single_cpm_fem_ode(dimval,mcsval,i,spreadval,tmaxval,Jhalf,dividing);
pause(30);
end
toc
%% size histo

figure
for i =10:10:1000
dat1=csize{1,i};
dat2=reshape
xlim([0 90])
title(['MCS: ' num2str(i)])
pause(.1)
end