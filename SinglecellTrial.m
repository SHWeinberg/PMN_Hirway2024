%Shreyas Hirway 6/9/19
%% Initial given parameters of cell in E state
close all
clc
init=struct;
einit=[0.16 0.01 0.01 0.38 0.03 0.01 0.35 3.2 0];
minit=[2.0696 0.3098 2.6581 0.0352 0.2772 2.7956 0.0075 0.0485 3.1515];
c=einit;
init.T0= c(1);% T TGF-B
init.s0=c(2);% s snail
init.S0=c(3);% S SNAIL
init.R30=c(4);% R3 miR-34
init.z0=c(5);% z zeb
init.Z0=c(6);% Z ZEB
init.R20=c(7);% R2 miR-200
init.E0=c(8);% E E-cadherin
init.N0=c(9);% N N-cadherin
init.exo=0; % Exogenous TGFB in the system
init.Tbas=0.06;
init.mir2=0.0002;
init.scale=1;
init.conv=1;

init.time=[0:1:24*20];% hours
t=init.time;
%%
init.exo=1.8;
x=SingleCell(init);

%% Plotting Concentration Graphs (9)
%set(gcf,'Position',[700 100 1000 800])

figure
subplot(5,1,1);
plot(t,x(:,2)) % snail
title('Concentration of snail')
ylim([0,0.4])

subplot(5,1,2);
plot(t,x(:,3)) % SNAIL
title('Concentration of SNAIL')
ylim([-0.25,4])

subplot(5,1,3);
plot(t,x(:,4)) % miR34
title('Concentration of miR34')
ylim([0,0.4])

subplot(5,1,4);
plot(t,x(:,5)) % zeb
title('Concentration of zeb')
ylim([0,0.4])

subplot(5,1,5);
plot(t,x(:,6)) % ZEB
title('Concentration of ZEB')
ylim([-0.5,4])

figure
subplot(4,1,1);
plot(t,x(:,7)) % miR200
title('Concentration of miR200')
ylim([0,0.5])

subplot(4,1,2);
plot(t,x(:,8)) % E Cadherin
title('Concentration of E cadherin')
ylim([-0.25,4])

subplot(4,1,3);
plot(t,x(:,9)) % N Cadherin
title('Concentration of N Cadherin')
ylim([-0.25,4])

subplot(4,1,4);
plot(t,x(:,1)) % TGF-B
title('Concentration of TGF-B')
xlabel('Time (hours)')
ylim([-.5,2.5])


