% For Testing Regular CPM with ECM stuff
dimval=100;
mcsval=150;
tgfbval=0;
spreadval=0.015*4;
tmaxval=0;
Jmin=4819*19;
dividing=1.25;%/4;
al=0.05;
n=2;
Jhalf0=Jmin/((al^-1-1)^(1/n));
Tmax=2*.125;
tmax0=1*.125;
Jmin= Jhalf0*(al^-1 -1) ^(1/n) / (2*Tmax/tmax0 - 1)^(1/n);
Jhalf=Jmin/((al^-1-1)^(1/n));
comment='';
pZEB=0.5;
pR200=1.5;
tmaxval=Tmax;
extrascale1=0.9;
extrascale2=5;

(9574.38-Jhalf)/9574.38;
intrazeb1=1;
intrasnail2=1;
kdae_scale=1;
%%
mcsval=10;
clc
run_single_cpm_fem_ode(dimval,mcsval,tgfbval,spreadval,tmaxval,Jhalf,dividing,extrascale1,extrascale2,intrazeb1, intrasnail2,pZEB,pR200,kdae_scale,comment);
pause(1)
close all
%%
% run j = 14 16 18, if it doesnt work, run 15,16- SUH010321
%%
pZEB=1;
pR200=1;
intrazeb1=1;
intrasnail2=1;

kdae_scale=1;
comment='JhalfRR';
%tgf_range=[0 5];
jhalf_range=[7 11 19 23];
for mcsval=1201:1205

for tv=[0.15] % and 0.15
for kdae_scale=1
    
    for extra1=[0.96]
        for extra2=[0.6]
            for j=[7 11 19 23]
             %parfor tgf_val1= 1:2
                tgf=5; % will only run 5 tgf_range(tgf_val1);
%                 j=jhalf_range(j_index);
                Jmin=4819*j;
                al=0.05;
                n=2;
                Jhalf0=Jmin/((al^-1-1)^(1/n));
                Tmax=2*tv;
                tmax0=1*tv;
                Jmin= Jhalf0*(al^-1 -1) ^(1/n) / (2*Tmax/tmax0 - 1)^(1/n);
                Jhalf=Jmin/((al^-1-1)^(1/n));
                tmaxval=Tmax;
                fac=[0.2];
                pZEB=1-fac;
                pR200=1+fac;
                extrascale1=extra1;
                extrascale2=extra2;
                tgfbval=tgf;
                disp(['Running tgfb: ' num2str(tgfbval)]);
                run_single_cpm_fem_ode(dimval,mcsval,tgfbval,spreadval,tmaxval,Jhalf,dividing,extrascale1,extrascale2,intrazeb1, intrasnail2,pZEB,pR200,kdae_scale,comment);
                pause(1)
                close all
            end
        end
    end
end
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