% Will allow for matrix of parameters, so that each cell can have its own
% set, creating heterogeneity- SUH 071921

function [init]=tian_model_baseline_params_v4(parammatrix,tmaxval,Jhalf,extrascale1,extrascale2,intrazeb1,intrasnail2,pZEB,pR200,kdae_scale,ke_scale,pop_index,NRc)
% parameters
% Compatible with CellFunc_v5
% New for F equation- SUH 112519

% New
init=struct;


%Extracellular concentrations defined here- SUH 111320
extra=struct;
extra.param_js=0.625;
extra.param_ke=0.525;
extra.param_ka=4.45;
extra.param_kdae=1.225;
 
init.scale=2; % OG was 2
init.conv=1; % OG was 1
 
 
%init.exo=3; % Exogenous TGFB in the system

 
% init.mask=mat;
% init.time=[0:1:24*20];% hours

 modelparameters=zeros(53,NRc);
%Scaling parameters for Partial state
partialscaleZEB=pZEB; % OG was 1, changed based on changingpartialstategraph
partialscalemir200=pR200; % OG was 1, changed based on changingpartialstategraph
%SUH-112320
 
% modelparameters(1,:)=0.06;  % Basal production rate of TGF-? 0.06 ?M/hr
% modelparameters(2,:)=1.2;  % Production rate of TGF-? 1.2 ?M/hr
% modelparameters(3,:)=0.06;  % Michaelis constant of miR200-dependent inhibition of TGF-? expression 0.06 ?M
% modelparameters(4,:)=0.6;  % Degradation rate of TGF-? 0.6/hr
%  
% modelparameters(5,:)=0.0006;  % Basal transcription rate of snail1 0.0006 ?M/hr
% modelparameters(6,:)=0.03;  %Transcription rate of snail1 0.03 ?M/hr
% modelparameters(7,:)=0.09*intrasnail2;  %Degradation rate of snail1 mRNA 0.09/hr
% modelparameters(8,:)=1.6;  %Michaelis constant of TGF-?-dependent snail1 translation 1.6 ?M
%  
% modelparameters(9,:)=17;  %Translation rate of snail1 mRNA 17 ?M/hr
% modelparameters(10,:)=1.66;  % Degradation rate of SNAIL1 1.66/hr
% modelparameters(11,:)=0.08;  % Michaelis constant of miR34-dependent inhibition of snail1 translation 0.08 ?M
%  
% modelparameters(12,:)=0.0012;  % Basal production rate of miR-34 0.0012 ?M/hr
% modelparameters(13,:)=0.012;  % Production rate of miR-34 0.012 ?M/hr
% modelparameters(14,:)=0.035;  % Degradation rate of miR-34 0.035/hr
% modelparameters(15,:)=0.15;  % Michaelis constant of SNAIL1-dependent inhibition of miR-34 production 0.15 ?M
% modelparameters(16,:)=0.36;  % Michaelis constant of ZEB-dependent inhibition of miR-34 production 0.36 ?M
%  
% modelparameters(17,:)=0.003;  % Basal transcription rate of zeb 0.003 ?M/hr
% modelparameters(18,:)=0.06;  % Transcription rate of zeb 0.06 ?M/hr
% modelparameters(19,:)=3.5;  % Michaelis constant of SNAIL1-dependent zeb transcription 3.5 ?M
% modelparameters(20,:)=0.09*intrazeb1;  % Degradation rate of zeb mRNA 0.09/hr
%  
% modelparameters(21,:)=1.66;  % Degradation rate of ZEB 1.66/hr
% modelparameters(22,:)=17*partialscaleZEB;  % Translation rate of zeb mRNA 17 ?M/hr *Could decrease this for partial state- SUH 111020
% modelparameters(23,:)=0.06;  % Michaelis constant of miR34-dependent inhibition of zeb mRNA translation 0.06 ?M
%  
% modelparameters(24,:)=0.0002;  % Basal production rate of miR-200 0.0002 ?M/hr
% modelparameters(25,:)=0.012*partialscalemir200;  % Production rate of miR-200 0.012 ?M/hr *Could increase this for partial state- SUH 111020
% modelparameters(26,:)=0.035;  % Degradation rate of miR-200 0.035/hr
% modelparameters(27,:)=5;  % Michaelis constant of SNAIL1-dependent inhibition of miR-200 production 5 ?M
% modelparameters(28,:)=0.2;  % Michaelis constant of ZEB-dependent inhibition of miR-200 production 0.2 ?M
%  
% modelparameters(29,:)=1;  % Production rate 1 of E-cadherin production 1 ?M/hr
% modelparameters(30,:)=0.6;  % Production rate 2 of E-cadherin production 0.6 ?M/hr
% modelparameters(31,:)=0.5;  % Degradation rate of E-cadherin 0.5/hr
% modelparameters(32,:)=0.2;  % Michaelis constant of SNAIL1-dependent inhibition of E-cadherin production 0.2 ?M
% modelparameters(33,:)=0.5;  % Michaelis constant of ZEB-dependent inhibition of E-cadherin production 0.5 ?M
%  
% modelparameters(34,:)=1;  % Production rate 1 of N-cadherin production 1 ?M/hr
% modelparameters(35,:)=0.6;  % Production rate 2 of N-cadherin production 0.6 ?M/hr
% modelparameters(36,:)=0.2;  % Michaelis constant of SNAIL1-dependent N-cadherin production 0.2?M
% modelparameters(37,:)=0.5;  % Michaelis constant of ZEB-dependent N-cadherin production 0.5 ?M
% modelparameters(38,:)=0.5;  % Degradation rate of N-cadherin 0.5/hr
%  
% modelparameters(39,:)=2;  % Hill coefficient of TGF-?-dependent SNAIL1 expression 2
% modelparameters(40,:)=2;  % Hill coefficient of miR-200-dependent inhibition 2
% modelparameters(41,:)=2;  % Hill coefficient of miR-34-dependent inhibition 2
% modelparameters(42,:)=2;  % Hill coefficient of ZEB-dependent inhibition 2
% modelparameters(43,:)=2;  % Hill coefficient of SNAIL1-dependent activation or inhibition 2
%  

for i=1:53
    modelparameters(i,:)=parammatrix(i,:);% set model parameters from parammatrix
end

% Special scaling factor cases
modelparameters(14,:)=modelparameters(14,:)*intrasnail2;
modelparameters(16,:)=modelparameters(16,:)*intrazeb1;
modelparameters(11,:)=modelparameters(11,:)*partialscaleZEB;
modelparameters(12,:)=modelparameters(12,:)*partialscalemir200;

modelparameters(50,:)=modelparameters(50,:)*extra.param_ka; %extra.ka2 is scaled by param_ka
%modelparameters(52,:)=modelparameters(52,:)*extra.param_kdae; % kdae2 scaled by param_kdae

init.parametermatrix=modelparameters;
%TGFB- added to parameter matrix in outer execution file- SUH 10/22/21
extra.kbT=10.656; % 0.008- calculated from 12th-14th paper
extra.kE1=0.006; % from 10-40nm ECM fibronectin prod
extra.kE2=0.0273;
%mdlparams.kE=((kE2-kE1).*Ncad_val/3.1515 + kE1)*init.param_ke;
extra.kdE=0.6; %from Tian Model
extra.fa=10; % Scaling Factor
extra.ka1=1/(extra.fa*12);
extra.ka2=(1/12)*extra.param_ka;
%mdlparams.ka=(ka2-ka1).*Ncad_val/3.1515 + ka1;
extra.fdAE=5; % Scaling Factor
extra.kdAE2=(720/500/25)*extra.param_kdae;
extra.kdAE1=extra.kdAE2/extra.fdAE; % from ECM degradation graph
%mdlparams.kdAE=(kdAE2-kdAE1).*Ncad_val/3.1515 + kdAE1;


extra.Jhalf=Jhalf;%1105.59636910808;
extra.nf=2;
extra.nf2=2;
extra.tmax=tmaxval; %Tmax=2.35, OG was fmax= Tmaxval/Js
extra.scaling1=extrascale1;
extra.scaling2=extrascale2;
% extra.intrazeb1=intrazeb1;
% extra.intrasnail2=intrasnail2;
extra.kdae_scale=kdae_scale;
 extra.ke_scale=ke_scale;
 
init.extra=extra;

end