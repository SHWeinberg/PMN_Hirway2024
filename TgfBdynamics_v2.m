%% Base Test Case
mat=zeros(10,10);
mat(4:7,6:7)=1;
mat(4:5,5)=1;
mat(7:9,1:2)=2;
mat(8,3)=2;
cells=max(max(mat));
cellparameters=ones(cells,1);
mat1=reshape(mat,size(mat,1)^2,1);
tgfbvals=[0:0.1:10];
endt=zeros(5,size(tgfbvals,2));
jsvals=zeros(5,1);
%%
diff1=0.001;
scale1=[1-2*diff1,1-diff1,1,1+diff1,1+diff1*2];
j=1;
tic
for j=1:5
for i=1:size(tgfbvals,2)
init=struct;
init.T0= 0.16;% T TGF-B        
init.s0=0.01;% s snail
init.S0=0.01;% S SNAIL
init.R30=0.38;% R3 miR-34
init.z0=0.03;% z zeb
init.Z0=0.01;% Z ZEB
init.R20=0.35;% R2 miR-200
init.E0=3.2;% E E-cadherin
init.N0=0;% N N-cadherin

%init.Tbas=0.06;
%init.mir2=0.0002;
init.scale=2;
init.conv=1;

init.solT=0.1571;
init.ecmT=0.0072;
init.solE=0.01;
init.asmE=2.4569*10^-5;

% % Trying to fix errors
% init.solT=0.1571;
% init.solE=0.0072;
% init.asmE=0.01;
% init.ecmT=2.4569*10^-5;

init.exo=tgfbvals(i); % Exogenous TGFB in the system


% Extracellular variables based on parameter study
init.param_js=0.625;%Tian0matched(i,2);
init.param_ke=0.525;%Tian0matched(i,3);
init.param_ka=4.45;%Tian0matched(i,4);
init.param_kdae=1.225;%Tian0matched(i,5);

init.mask=mat;
init.time=[0:1:24*20];% hours

 parametermatrix=zeros(43,size(cellparameters,1));
%Scaling parameters for Partial state
partialscaleZEB=1;
partialscalemir200=1;


parametermatrix(1,:)=0.06;  % Basal production rate of TGF-? 0.06 ?M/hr
parametermatrix(2,:)=1.2;  % Production rate of TGF-? 1.2 ?M/hr
parametermatrix(3,:)=0.06;  % Michaelis constant of miR200-dependent inhibition of TGF-? expression 0.06 ?M
parametermatrix(4,:)=0.6;  % Degradation rate of TGF-? 0.6/hr
 
parametermatrix(5,:)=0.0006;  % Basal transcription rate of snail1 0.0006 ?M/hr
parametermatrix(6,:)=0.03;  %Transcription rate of snail1 0.03 ?M/hr
parametermatrix(7,:)=0.09*scale1(j);  %Degradation rate of snail1 mRNA 0.09/hr
parametermatrix(8,:)=1.6;  %Michaelis constant of TGF-?-dependent snail1 translation 1.6 ?M
 
parametermatrix(9,:)=17;  %Translation rate of snail1 mRNA 17 ?M/hr
parametermatrix(10,:)=1.66;  % Degradation rate of SNAIL1 1.66/hr
parametermatrix(11,:)=0.08;  % Michaelis constant of miR34-dependent inhibition of snail1 translation 0.08 ?M
 
parametermatrix(12,:)=0.0012;  % Basal production rate of miR-34 0.0012 ?M/hr
parametermatrix(13,:)=0.012;  % Production rate of miR-34 0.012 ?M/hr
parametermatrix(14,:)=0.035;  % Degradation rate of miR-34 0.035/hr
parametermatrix(15,:)=0.15;  % Michaelis constant of SNAIL1-dependent inhibition of miR-34 production 0.15 ?M
parametermatrix(16,:)=0.36;  % Michaelis constant of ZEB-dependent inhibition of miR-34 production 0.36 ?M
 
parametermatrix(17,:)=0.003;  % Basal transcription rate of zeb 0.003 ?M/hr
parametermatrix(18,:)=0.06;  % Transcription rate of zeb 0.06 ?M/hr
parametermatrix(19,:)=3.5;  % Michaelis constant of SNAIL1-dependent zeb transcription 3.5 ?M
parametermatrix(20,:)=0.09;  % Degradation rate of zeb mRNA 0.09/hr
 
parametermatrix(21,:)=1.66;  % Degradation rate of ZEB 1.66/hr
parametermatrix(22,:)=17*partialscaleZEB;  % Translation rate of zeb mRNA 17 ?M/hr *Could decrease this for partial state- SUH 111020
parametermatrix(23,:)=0.06;  % Michaelis constant of miR34-dependent inhibition of zeb mRNA translation 0.06 ?M
 
parametermatrix(24,:)=0.0002;  % Basal production rate of miR-200 0.0002 ?M/hr
parametermatrix(25,:)=0.012*partialscalemir200;  % Production rate of miR-200 0.012 ?M/hr *Could increase this for partial state- SUH 111020
parametermatrix(26,:)=0.035;  % Degradation rate of miR-200 0.035/hr
parametermatrix(27,:)=5;  % Michaelis constant of SNAIL1-dependent inhibition of miR-200 production 5 ?M
parametermatrix(28,:)=0.2;  % Michaelis constant of ZEB-dependent inhibition of miR-200 production 0.2 ?M
 
parametermatrix(29,:)=1;  % Production rate 1 of E-cadherin production 1 ?M/hr
parametermatrix(30,:)=0.6;  % Production rate 2 of E-cadherin production 0.6 ?M/hr
parametermatrix(31,:)=0.5;  % Degradation rate of E-cadherin 0.5/hr
parametermatrix(32,:)=0.2;  % Michaelis constant of SNAIL1-dependent inhibition of E-cadherin production 0.2 ?M
parametermatrix(33,:)=0.5;  % Michaelis constant of ZEB-dependent inhibition of E-cadherin production 0.5 ?M
 
parametermatrix(34,:)=1;  % Production rate 1 of N-cadherin production 1 ?M/hr
parametermatrix(35,:)=0.6;  % Production rate 2 of N-cadherin production 0.6 ?M/hr
parametermatrix(36,:)=0.2;  % Michaelis constant of SNAIL1-dependent N-cadherin production 0.2?M
parametermatrix(37,:)=0.5;  % Michaelis constant of ZEB-dependent N-cadherin production 0.5 ?M
parametermatrix(38,:)=0.5;  % Degradation rate of N-cadherin 0.5/hr
 
parametermatrix(39,:)=2;  % Hill coefficient of TGF-?-dependent SNAIL1 expression 2
parametermatrix(40,:)=2;  % Hill coefficient of miR-200-dependent inhibition 2
parametermatrix(41,:)=2;  % Hill coefficient of miR-34-dependent inhibition 2
parametermatrix(42,:)=2;  % Hill coefficient of ZEB-dependent inhibition 2
parametermatrix(43,:)=2;  % Hill coefficient of SNAIL1-dependent activation or inhibition 2
 
%TGFB
extra=struct;
extra.kbT=10.656; % 0.008- calculated from 12th-14th paper
extra.kE1=0.006; % from 10-40nm ECM fibronectin prod
extra.kE2=0.0273;
%mdlparams.kE=((kE2-kE1).*Ncad_val/3.1515 + kE1)*init.param_ke;
extra.kdE=0.6; %from Tian Model
extra.fa=10; % Scaling Factor
extra.ka1=1/(extra.fa*12);
extra.ka2=(1/12)*init.param_ka;
%mdlparams.ka=(ka2-ka1).*Ncad_val/3.1515 + ka1;
extra.fdAE=5; % Scaling Factor
extra.kdAE2=(720/500/25)*init.param_kdae;
extra.kdAE1=extra.kdAE2/extra.fdAE; % from ECM degradation graph
%mdlparams.kdAE=(kdAE2-kdAE1).*Ncad_val/3.1515 + kdAE1;
init.extra=extra;
init.parameters=parametermatrix;

disp(['Running Cell Simulation: ' num2str(tgfbvals(i))])
[f1,extra1]=Multicell_v3(init,cells);
init.time=init.time/24;
f1cells=f1(:,301:316);

endt(j,i)=f1cells(end,end);
jsvals(j)=parametermatrix(6,1);
end
end
toc

%
figure
plot(tgfbvals, endt(1,:));
 hold on
 plot(tgfbvals, endt(2,:));
 hold on
 plot(tgfbvals, endt(3,:));
hold on
plot(tgfbvals, endt(4,:));
hold on
plot(tgfbvals, endt(5,:));
xlabel('TGFB');
ylabel('N-cad');
legend(num2str(scale1(1)),num2str(scale1(2)),num2str(scale1(3)),num2str(scale1(4)),num2str(scale1(5)));
%title(['Parameter : ' num2str(parametermatrix(2,1))]);
findfigs
beep
%%
figure
plot(init.time,f1cells(:,15));
hold on
plot(init.time,f1cells(:,16));
title('Cell Ncad vs Time');
xlabel('Time');
ylabel('Ncad');

figure
solT=extra1(:,1:100);
asmE=extra1(:,101:200);
ecmT=extra1(:,201:300);
subplot(1,3,1)
plot(init.time,solT(:,8));
xlabel('Time');
ylabel('solT');
subplot(1,3,2)
plot(init.time,asmE(:,8));
xlabel('Time');
ylabel('asmE');
subplot(1,3,3)
plot(init.time,ecmT(:,8));
xlabel('Time');
ylabel('ecmT');


findfigs






%% ************************************************************************
figure;
axes('NextPlot','add');
%% Ncad vs Tgfb Graph
exovals=[0:0.05:4];
Ncadvals=zeros(size(exovals));
for i=1:size(exovals,2)
    tic
mat=zeros(10,10);
mat(4:7,6:7)=1;
mat(4:5,5)=1;
mat(7:9,1:2)=2;
mat(8,3)=2;
cells=max(max(mat));
cellparameters=ones(cells,1);
mat1=reshape(mat,size(mat,1)^2,1);
init=struct;
init.T0= 0.16;% T TGF-B        
init.s0=0.01;% s snail
init.S0=0.01;% S SNAIL
init.R30=0.38;% R3 miR-34
init.z0=0.03;% z zeb
init.Z0=0.01;% Z ZEB
init.R20=0.35;% R2 miR-200
init.E0=3.2;% E E-cadherin
init.N0=0;% N N-cadherin

init.Tbas=0.06;
init.mir2=0.0002;
init.scale=2;
init.conv=1;

init.solT=0.1571;
init.ecmT=0.0072;
init.solE=0.01;
init.asmE=2.4569*10^-5;

init.exo=exovals(i); % Exogenous TGFB in the system
% Extracellular variables based on parameter study
init.param_js=0.625;%Tian0matched(i,2);
init.param_ke=0.525;%Tian0matched(i,3);
init.param_ka=4.45;%Tian0matched(i,4);
init.param_kdae=1.225;%Tian0matched(i,5);

init.mask=mat;
init.time=[0:1:24*20];% hours


mdlparams=struct;

%Scaling parameters for Partial state
partialscaleZEB=1;
partialscalemir200=1;

mdlparams.k0T=0.06*cellparameters;  % Basal production rate of TGF-? 0.06 ?M/hr
mdlparams.kT=1.2*cellparameters;  % Production rate of TGF-? 1.2 ?M/hr
mdlparams.JT=0.06*cellparameters;  % Michaelis constant of miR200-dependent inhibition of TGF-? expression 0.06 ?M
mdlparams.kdT=0.6*cellparameters;  % Degradation rate of TGF-? 0.6/hr

mdlparams.k0s=0.0006*cellparameters;  % Basal transcription rate of snail1 0.0006 ?M/hr
mdlparams.ks=0.03*cellparameters;  %Transcription rate of snail1 0.03 ?M/hr
mdlparams.kds=0.09*cellparameters;  %Degradation rate of snail1 mRNA 0.09/hr
mdlparams.Js=1.6*cellparameters;  %Michaelis constant of TGF-?-dependent snail1 translation 1.6 ?M

mdlparams.kS=17*cellparameters;  %Translation rate of snail1 mRNA 17 ?M/hr
mdlparams.kdS=1.66*cellparameters;  % Degradation rate of SNAIL1 1.66/hr
mdlparams.JS=0.08*cellparameters;  % Michaelis constant of miR34-dependent inhibition of snail1 translation 0.08 ?M

mdlparams.k03=0.0012*cellparameters;  % Basal production rate of miR-34 0.0012 ?M/hr
mdlparams.k3=0.012*cellparameters;  % Production rate of miR-34 0.012 ?M/hr
mdlparams.kd3=0.035*cellparameters;  % Degradation rate of miR-34 0.035/hr
mdlparams.J13=0.15*cellparameters;  % Michaelis constant of SNAIL1-dependent inhibition of miR-34 production 0.15 ?M
mdlparams.J23=0.36*cellparameters;  % Michaelis constant of ZEB-dependent inhibition of miR-34 production 0.36 ?M

mdlparams.k0z=0.003*cellparameters;  % Basal transcription rate of zeb 0.003 ?M/hr
mdlparams.kz=0.06*cellparameters;  % Transcription rate of zeb 0.06 ?M/hr
mdlparams.Jz=3.5*cellparameters;  % Michaelis constant of SNAIL1-dependent zeb transcription 3.5 ?M
mdlparams.kdz=0.09*cellparameters;  % Degradation rate of zeb mRNA 0.09/hr

mdlparams.kdZ=1.66*cellparameters;  % Degradation rate of ZEB 1.66/hr
mdlparams.kZ=17*partialscaleZEB*cellparameters;  % Translation rate of zeb mRNA 17 ?M/hr *Could decrease this for partial state- SUH 111020
mdlparams.JZ=0.06*cellparameters;  % Michaelis constant of miR34-dependent inhibition of zeb mRNA translation 0.06 ?M

mdlparams.k02=0.0002*cellparameters;  % Basal production rate of miR-200 0.0002 ?M/hr
mdlparams.k2=0.012*partialscalemir200*cellparameters;  % Production rate of miR-200 0.012 ?M/hr *Could increase this for partial state- SUH 111020
mdlparams.kd2=0.035*cellparameters;  % Degradation rate of miR-200 0.035/hr
mdlparams.J12=5*cellparameters;  % Michaelis constant of SNAIL1-dependent inhibition of miR-200 production 5 ?M
mdlparams.J22=0.2*cellparameters;  % Michaelis constant of ZEB-dependent inhibition of miR-200 production 0.2 ?M

mdlparams.ke1=1*cellparameters;  % Production rate 1 of E-cadherin production 1 ?M/hr
mdlparams.ke2=0.6*cellparameters;  % Production rate 2 of E-cadherin production 0.6 ?M/hr
mdlparams.kde=0.5*cellparameters;  % Degradation rate of E-cadherin 0.5/hr
mdlparams.J1e=0.2*cellparameters;  % Michaelis constant of SNAIL1-dependent inhibition of E-cadherin production 0.2 ?M
mdlparams.J2e=0.5*cellparameters;  % Michaelis constant of ZEB-dependent inhibition of E-cadherin production 0.5 ?M

mdlparams.kn1=1*cellparameters;  % Production rate 1 of N-cadherin production 1 ?M/hr
mdlparams.kn2=0.6*cellparameters;  % Production rate 2 of N-cadherin production 0.6 ?M/hr
mdlparams.J1n=0.2*cellparameters;  % Michaelis constant of SNAIL1-dependent N-cadherin production 0.2?M
mdlparams.J2n=0.5*cellparameters;  % Michaelis constant of ZEB-dependent N-cadherin production 0.5 ?M
mdlparams.kdn=0.5*cellparameters;  % Degradation rate of N-cadherin 0.5/hr

mdlparams.nt=2*cellparameters;  % Hill coefficient of TGF-?-dependent SNAIL1 expression 2
mdlparams.nr2=2*cellparameters;  % Hill coefficient of miR-200-dependent inhibition 2
mdlparams.nr3=2*cellparameters;  % Hill coefficient of miR-34-dependent inhibition 2
mdlparams.nz=2*cellparameters;  % Hill coefficient of ZEB-dependent inhibition 2
mdlparams.ns=2*cellparameters;  % Hill coefficient of SNAIL1-dependent activation or inhibition 2

%TGFB
mdlparams.kbT=10.656; % 0.008- calculated from 12th-14th paper
mdlparams.kE1=0.006; % from 10-40nm ECM fibronectin prod
mdlparams.kE2=0.0273;
%mdlparams.kE=((kE2-kE1).*Ncad_val/3.1515 + kE1)*init.param_ke;
mdlparams.kdE=0.6; %from Tian Model
mdlparams.fa=10; % Scaling Factor
mdlparams.ka1=1/(mdlparams.fa*12);
mdlparams.ka2=(1/12)*init.param_ka;
%mdlparams.ka=(ka2-ka1).*Ncad_val/3.1515 + ka1;
mdlparams.fdAE=5; % Scaling Factor
mdlparams.kdAE2=(720/500/25)*init.param_kdae;
mdlparams.kdAE1=mdlparams.kdAE2/mdlparams.fdAE; % from ECM degradation graph
%mdlparams.kdAE=(kdAE2-kdAE1).*Ncad_val/3.1515 + kdAE1;

init.parameters=mdlparams;

%disp(['Running Cell Simulation'])
f1=Multicell_v2(init,cells);
init.time=init.time/24;
f1cells=f1(:,401:416);


% Finding comparisons to Tian Model 
%f1cells(:,15)- Ncad
% Day 10- 240, Day 15- 360

% if ((f1cells(240,15)<2.5) && (f1cells(360,15)>3.14))
%     disp(['Comparable conditions: ' num2str(init.exo) ' ExoT'])
%     disp(['and ' num2str(init.param_js) ' Js'])
%     disp(['and ' num2str(init.param_ke) ' Ke'])
%     disp(['and ' num2str(init.param_ka) ' Ka'])
%     disp(['and ' num2str(init.param_kdae) ' Kdae'])
%     Conditionmatched(count,1)=init.exo;
%     Conditionmatched(count,2)=init.param_js;
%     Conditionmatched(count,3)=init.param_ke;
%     Conditionmatched(count,4)=init.param_ka;
%     Conditionmatched(count,5)=init.param_kdae;
%     count=count+1;
%     
% else 
%     disp('Conditions didnt match')
% end

Ncadvals(i)=f1cells(end,16);
end

toc
plot(exovals,Ncadvals, 'DisplayName', ['ZEB: ' num2str(partialscaleZEB) ', miR200: ' num2str(partialscalemir200)]);
%title();
%%
legend('show');
%% Extracellular Concentration Graphs
f1mat=f1(:,1:400);
figure
subplot(2,2,1)
plot(init.time,f1mat(:,7))
ylabel('solT');
subplot(2,2,2)
plot(init.time,f1mat(:,107))
ylabel('solE');
subplot(2,2,3)
plot(init.time,f1mat(:,207))
ylabel('asmE');
subplot(2,2,4)
plot(init.time,f1mat(:,307))
ylabel('ecmT');

% Visualizing Cell Dynamics

% figure
% subplot(4,2,1)
% plot(init.time,f1cells(:,1));
% ylabel('snail1');
% ylim([0,0.4]);
% subplot(4,2,2)
% plot(init.time,f1cells(:,2));
% ylim([0,0.4]);
% subplot(4,2,3)
% plot(init.time,f1cells(:,3));
% ylabel('SNAIL1');
% ylim([0,2]);
% subplot(4,2,4)
% plot(init.time,f1cells(:,4));
% ylim([0,2]);
% subplot(4,2,5)
% plot(init.time,f1cells(:,5));
% ylabel('mir34');
% ylim([0,0.4]);
% subplot(4,2,6)
% plot(init.time,f1cells(:,6));
% ylim([0,0.4]);
% subplot(4,2,7)
% plot(init.time,f1cells(:,7));
% ylabel('zeb');
% ylim([0,0.4]);
% subplot(4,2,8)
% plot(init.time,f1cells(:,8));
% ylim([0,0.4]);
% title(['Run: ' num2str(i)]);
% figure
% subplot(4,2,1)
% plot(init.time,f1cells(:,9));
% ylabel('ZEB');
% ylim([0,4]);
% subplot(4,2,2)
% plot(init.time,f1cells(:,10));
% ylim([0,4]);
% subplot(4,2,3)
% plot(init.time,f1cells(:,11));
% ylabel('miR200');
% ylim([0,0.4]);
% subplot(4,2,4)
% plot(init.time,f1cells(:,12));
% ylim([0,0.4]);
% subplot(4,2,5)
% plot(init.time,f1cells(:,13));
% ylabel('E-Cad');
% ylim([0,4]);
% subplot(4,2,6)
% plot(init.time,f1cells(:,14));
% ylim([0,4]);
% subplot(4,2,7)
% plot(init.time,f1cells(:,15));
% ylabel('N-Cad');
% ylim([0,4]);
% subplot(4,2,8)
% plot(init.time,f1cells(:,16));
% ylim([0,4]);
% title(['Run: ' num2str(i)]);

%% Visualizing Grid with Extracellular dynamics
f1mat=f1(:,1:400);
figure
for i=1:1:481
    %disp([num2str(i)]);
    maxv1=max(max(f1mat(:,1:100)));
    maxv2=max(max(f1mat(:,101:200)));
    maxv3=max(max(f1mat(:,201:300)));
    maxv4=max(max(f1mat(:,301:400)));
    
    matrx1=f1mat(i,1:100);
    matrx1=reshape(matrx1,10,10);
    matrx2=f1mat(i,101:200);
    matrx2=reshape(matrx2,10,10);
    matrx3=f1mat(i,201:300);
    matrx3=reshape(matrx3,10,10);
    matrx4=f1mat(i,301:400);
    matrx4=reshape(matrx4,10,10);
    
    subplot(2,2,1)
    imagesc(matrx1);   
    caxis([0 maxv1])
    
    subplot(2,2,2)
    imagesc(matrx2);   
    caxis([0 maxv2])
    
    subplot(2,2,3)
    imagesc(matrx3);   
    caxis([0 maxv3])
    
    subplot(2,2,4)
    imagesc(matrx4);   
    caxis([0 maxv4])
    pause(0.5)
    
end