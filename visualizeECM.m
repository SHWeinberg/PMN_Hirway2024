%% ECM plotting
%% Visualize ECM Snapshots 7/15 *******************************************- Used in CPM paper for ECM visualization
tic
count1=17;
scalev=1.35;
for paramval=[0.1]
filename=['test_100p_1402mcs_1tgfb__Random_divf_0.5_param_14_scale_0.5ECMscaled' num2str(paramval) '_S&K'];

load([filename '.mat']);
%load('test_100p_2502mcs_1tgfb_0.075spread1.25div_extraA0.96_extraBT0.6_ZEB1.125_SNAIL1_kDAE.scale1_base.mat');
%
ReadObj = VideoReader(['Cell_Phenotype_' filename '.mp4']);
ReadObj2 = VideoReader(['Cell_Types_' filename '.mp4']);

grid=100;
fulltime=size(csize,2);
solTvar=zeros(grid*grid,fulltime);
solEvar=zeros(grid*grid,fulltime);
asmEvar=zeros(grid*grid,fulltime);
ecmTvar=zeros(grid*grid,fulltime);
for MCS=1:fulltime
    
    maskG=ctag(:,MCS);
    stateV=statevars{MCS};
    NRC=size(csize{MCS},1);
    k0T_val=zeros(1,NRC);
    k0T_val(:)=0.06;
    kT_val=zeros(1,NRC);
    kT_val(:)=1.2;
    nx=grid*grid;
    MR200_val=zeros(sqrt(nx),sqrt(nx));
    R2Val=stateV((grid*grid*3)+(5*NRC)+1: (grid*grid*3)+(6*NRC));
    for i=1:NRC
        M200_temp=R2Val(i);
        MR200_val(find(maskG==i))=M200_temp;
    end
    MR200_val=reshape(MR200_val,nx,1);
    
    JT_val=zeros(1,NRC);
    JT_val(:)=0.06;
    nr2_val=zeros(1,NRC);
    nr2_val(:)=2;
    kT200 = k0T_val + kT_val./(1 + ([MR200_val]./JT_val).^nr2_val);
    kT200(isnan(kT200))=0;
    
    solEv=stateV(1:(grid*grid));
    totalEv=stateV((grid*grid+1):(grid*grid*2));
    totalTv=stateV((grid*grid*2+1):(grid*grid*3));
    kbT_mask=zeros(sqrt(nx),sqrt(nx));
    kbT_mask(:)=10.656*0.6; %par.extra.kbT*extrascaling2;
    Ncad_val=zeros(sqrt(nx),sqrt(nx));
    NVal=stateV((grid*grid*3)+(7*NRC)+1: (grid*grid*3)+(8*NRC));
    for i=1:NRC
        Ncad_temp=NVal(i);
        Ncad_val(find(maskG==i))=Ncad_temp;
    end
    Ncad_val=reshape(Ncad_val,nx,1);
    
    param_kdae=1.225;
    
    fdAE=5; % Scaling Factor
    extra.kdAE2=(720/500/25)*param_kdae;
    extra.kdAE1=extra.kdAE2/fdAE; % from ECM degradation graph
    %kdAE2=(720/500/25)*param_kdae;
    %kdAE1=kdAE2/fdAE; % from ECM degradation graph
    %kdAE=(kdAE2-kdAE1).*Ncad_val/3.1515 + kdAE1; % Uses Ncad so calculated here
    %kdAE=kdAE*1;
    kdAE2=parammatrix{MCS}(52,:)';%(720/500/25)*param_kdae;
    kdAE1=parammatrix{MCS}(53,:)';%kdAE2/fdAE; % from ECM degradation graph
    kdAE1_mask=zeros(sqrt(nx),sqrt(nx));
    kdAE1_mask(:)=extra.kdAE1;
    kdAE2_mask=zeros(sqrt(nx),sqrt(nx));
    
    for i=1:NRC
        %     kbT_temp=kbT(i);
        %     kbT_mask(find(maskG==i))=kbT_temp;
        
        kdAE1_temp=kdAE1(i);
        kdAE1_mask(find(maskG==i))=kdAE1_temp;
        
        kdAE2_temp=kdAE2(i);
        kdAE2_mask(find(maskG==i))=kdAE2_temp;
    end
    kdAE1_mask=reshape(kdAE1_mask,nx,1);
    kdAE2_mask=reshape(kdAE2_mask,nx,1);
    kbT_mask=reshape(kbT_mask,nx,1);
    kbT=kbT_mask;
    kdAE=(kdAE2_mask-kdAE1_mask).*Ncad_val/3.1515 + kdAE1_mask; % Uses Ncad so calculated here
    kdae_scale=1;
    kdAE=kdAE*kdae_scale;
    
    %ecmTv=(kbT.*(totalTv+totalEv)+kdAE - sqrt( (kbT.*(totalTv+totalEv)+kdAE).^2 - 4.*kbT.^2.*totalTv.*totalEv))./(2*kbT);
    ecmTv=(kbT.*(totalTv+totalEv)+kdAE.*(sign(maskG)) - sqrt( (kbT.*(totalTv+totalEv)+kdAE.*(sign(maskG))).^2 - 4.*kbT.^2.*totalTv.*totalEv))./(2*kbT); % added .*(sign(maskG)) to kdAE- SUH 10/25/21
    
    solTv=totalTv-ecmTv;
    asmEv=totalEv-ecmTv;
    
    % storing variables
    solTvar(:,MCS)=solTv;
    solEvar(:,MCS)=solEv;
    asmEvar(:,MCS)=asmEv;
    ecmTvar(:,MCS)=ecmTv;
    
end

%
cellloc=maskG>0;
totalvals=[solTvar(:,end),solEvar(:,end),asmEvar(:,end),ecmTvar(:,end)];
midvals=totalvals(1225,:);
avgvals=mean(totalvals(cellloc,:));

ECMrange(count1,1)=min(min(ecmTvar));
ECMrange(count1,2)=max(max(ecmTvar));
ECMrange(count1,3)=scalev;

toc
% 90% EcmT
% max90pct=maxEcmTvalue*0.9;
% avgT=zeros(fulltime,1);
% max90T=zeros(fulltime,1);
% for i=1:fulltime
%     temp=ecmTvar(:,i)>=max90pct;
%     max90T(i)=sum(sum(temp));
%     avgT(i)=mean(mean(ecmTvar(:,i)));
% end
% midvals

%Creating cell type mask
% celltypemask=zeros(size(ctag,1),size(ctag,2));
% for mcs=1:size(ctag,2)
%     for c=1:size(csize{mcs},1)
%         celltypemask(find(ctag(:,mcs)==c),mcs)=pop_index{mcs}(c);
%     end
% end

%
numFrames = ceil(ReadObj.FrameRate*ReadObj.Duration);
f=figure;
colormap(jet);
ReadObj.CurrentTime=0;
ReadObj2.CurrentTime=0;
MCSval=1;
vidfile = VideoWriter([filename 'cells&EcmT&type.mp4'],'MPEG-4');
open(vidfile);
while hasFrame(ReadObj)
    vidFrame = readFrame(ReadObj);
    subplot(1,3,1)
    imshow(vidFrame)
    
    subplot(1,3,2)
    imagesc(rot90(reshape(ecmTvar(:,MCSval),grid,grid)));
    set(gca,'YColor','white');
    set(gca,'XColor','white');
    title(['EcmT at MCS: ' num2str(MCSval)]);
    
    caxis([min(min(ecmTvar)), max(max(ecmTvar))]);
    colorbar
    
    vidFrame2 = readFrame(ReadObj2);
    subplot(1,3,3)
    imshow(vidFrame2)
    
    
    
    MCSval=MCSval+1;
    
    f.Position = [2,150,1916,459];
    F = getframe(f);
    [RGB] = frame2im(F);
    writeVideo(vidfile, RGB);
    pause(1/ReadObj.FrameRate);
end
close(vidfile)
close all

toc
end
%% Store in variables
tic
PMNrandom{1,1}='0.1tgfb';
PMNrandom{1,2}=avgT;
PMNrandom{1,3}=max90T;
toc
%% 1x3 for 90%, average and final EcmT
load('PMNtest1.mat')
figure
subplot(1,3,1)
plot(PMNtest1{1,2}, 'LineWidth', 2);
hold on
plot(PMNtest1{2,2},'LineWidth',2);
hold on
plot(PMNtest1{3,2}, 'LineWidth', 2);
hold on
plot(PMNtest1{4,2}, 'LineWidth', 2);

xlabel('Time (MCS)','FontSize',24);
ylabel('ECM-Bound TGF-?  (?M)','FontSize',24);
title('Average ECM-Bound TGF-?','FontSize',26);
ax = gca;
ax.FontSize = 18;

subplot(1,3,2)
plot((PMNtest1{1,3}/(100*100))*100, 'LineWidth', 2);
hold on
plot((PMNtest1{2,3}/(100*100))*100, 'LineWidth', 2);
hold on
plot((PMNtest1{3,3}/(100*100))*100, 'LineWidth', 2);
hold on
plot((PMNtest1{4,3}/(100*100))*100, 'LineWidth', 2);
legend(PMNtest1{1,1},PMNtest1{2,1},PMNtest1{3,1},PMNtest1{4,1},'FontSize',18);
xlabel('Time (MCS)','FontSize',24);
ylabel('Grid Percentage (%)','FontSize',24);
title('90% Max ECM-Bound TGF-?','FontSize',26);
ax = gca;
ax.FontSize = 18;

subplot(1,3,3)
plot([0, 0.1, 0.5, 1],[PMNtest1{1,2}(end),PMNtest1{2,2}(end),PMNtest1{3,2}(end),PMNtest1{4,2}(end)], 'LineWidth', 2);
xlabel('TGF-? Dose (?M)','FontSize',24);
ylabel('ECM-Bound TGF-?  (?M)','FontSize',24);
title('Final ECM-Bound TGF-?','FontSize',26);
ax = gca;
ax.FontSize = 18;
%% Plotting at certain MCS for 4 concentrations
figure
colormap(jet);
MCSrange=[];
for i=1:5
    MCSval=510+(i-1)*50;
    subplot(4,5,i)
    imagesc(rot90(reshape(solTvar(:,MCSval),grid,grid)));
    caxis([min(min(solTvar)), max(max(solTvar))]);
    colorbar
    title(['solT at ' num2str(MCSval)]);
    subplot(4,5,i+5)
    imagesc(rot90(reshape(solEvar(:,MCSval),grid,grid)));
    caxis([min(min(solEvar)), max(max(solEvar))]);
    colorbar
    title(['solE at' num2str(MCSval)]);
    subplot(4,5,i+10)
    imagesc(rot90(reshape(asmEvar(:,MCSval),grid,grid)));
    caxis([min(min(asmEvar)), max(max(asmEvar))]);
    colorbar
    title(['asmE at' num2str(MCSval)]);
    subplot(4,5,i+15)
    imagesc(rot90(reshape(ecmTvar(:,MCSval),grid,grid)));
    caxis([min(min(ecmTvar)), max(max(ecmTvar))]);
    colorbar
    title(['ecmT at' num2str(MCSval)]);
end
sgtitle('5 TGFB');
%% Plotting ecmT

figure
colormap(jet);
totalEvar=asmEvar+ecmTvar;
MCSrange=[1000 1500 2000]; % 1 figure with 600,700,800,900 and 1 with 1000,1100,1200,1300
for i=1:size(MCSrange,2)
    
    MCSval=MCSrange(i);
    
    subplot(1,size(MCSrange,2),i)
    imagesc(rot90(reshape(ecmTvar(:,MCSval),grid,grid)));
    caxis([min(min(ecmTvar)), max(max(ecmTvar))]);
    
    
    title(['ecmT at ' num2str(MCSval*4.8/60) ' Hrs']);
    
end
colorbar


%% save ecmt video alongside
tic
numFrames = ceil(ReadObj.FrameRate*ReadObj.Duration);
f=figure;
colormap(jet);
ReadObj.CurrentTime=0;
MCSval=1;
vidfile = VideoWriter([filename 'redo' 'cells&EcmT.mp4'],'MPEG-4');
open(vidfile);
while hasFrame(ReadObj)
    vidFrame = readFrame(ReadObj);
    subplot(1,3,1)
    imshow(vidFrame)
    
    subplot(1,3,2)
    imagesc(rot90(reshape(ecmTvar(:,MCSval),grid,grid)));
    set(gca,'YColor','white');
    set(gca,'XColor','white');
    title(['EcmT at MCS: ' num2str(MCSval)]);
    
    caxis([min(min(ecmTvar)), max(max(ecmTvar))]);
    colorbar
    subplot(1,3,3)
    imagesc(rot90(reshape(celltypemask(:,MCSval),grid,grid)));
    set(gca,'YColor','white');
    set(gca,'XColor','white');
    title(['Cell Types at MCS: ' num2str(MCSval)]);
    hAxes = gca;
    colour0 = [1 1 1];
    colour1 = [0 0.4470 0.7410];
    colour2= [0.6350 0.0780 0.1840];
    colormap( hAxes, [colour0; colour1; colour2] )
    
    colorbar('Ticks',[0,1,2],...
        'TickLabels',{'Grid','Kupffer','Stellate'})
    
    MCSval=MCSval+1;
    
    f.Position = [2,150,1916,459];
    F = getframe(f);
    [RGB] = frame2im(F);
    writeVideo(vidfile, RGB);
    pause(1/ReadObj.FrameRate);
end
close(vidfile)
toc