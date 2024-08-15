
%%
tic
runid=1000284;
SimID=runid;

before=['test_100p_800mcs_5tgfb__BeforeTumor_divf_0.5_param_14_scale_0.5ECMscaled0_' num2str(SimID)];
after=['test_100p_1000mcs_0tgfb__AfterTumor_divf_0.75_param_14_scale_0.5ECMscaled0_' num2str(SimID)];
load([before '.mat']);
videocount=videocount+1;
ReadObj = VideoReader(['Cell_Phenotype_' before '.mp4']);
ReadObj2 = VideoReader(['Cell_Types_' before '.mp4']);

grid=100;
fulltime=size(csize,2);
beforetime=fulltime;
beforevector=cell(fulltime,3);
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

beforecelltype=pop_index;
beforecellpheno=phenotype;

cellloc=maskG>0;
totalvals=[solTvar(:,end),solEvar(:,end),asmEvar(:,end),ecmTvar(:,end)];
midvals=totalvals(1225,:);
avgvals=mean(totalvals(cellloc,:));

ReadObj.CurrentTime=0;
ReadObj2.CurrentTime=0;
MCSval=1;
beforeEcmT=[min(min(ecmTvar)), max(max(ecmTvar))];
while hasFrame(ReadObj)
    vidFrame = readFrame(ReadObj);
    beforevector{MCSval,1}=vidFrame;
    
    beforevector{MCSval,2}=rot90(reshape(ecmTvar(:,MCSval),grid,grid));
    
    vidFrame2 = readFrame(ReadObj2);
    beforevector{MCSval,3}=vidFrame2;
    MCSval=MCSval+1;
end

toc
% After tumor is added

load([after '.mat']);

%load('test_100p_2502mcs_1tgfb_0.075spread1.25div_extraA0.96_extraBT0.6_ZEB1.125_SNAIL1_kDAE.scale1_base.mat');
%
ReadObj = VideoReader(['Cell_Phenotype_' after '.mp4']);
ReadObj2 = VideoReader(['Cell_Types_' after '.mp4']);

grid=100;
fulltime=size(csize,2);
aftertime=fulltime;
aftervector=cell(fulltime,3);
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

aftercelltype=pop_index;
aftercellpheno=phenotype;

cellloc=maskG>0;
totalvals=[solTvar(:,end),solEvar(:,end),asmEvar(:,end),ecmTvar(:,end)];
midvals=totalvals(1225,:);
avgvals=mean(totalvals(cellloc,:));

ReadObj.CurrentTime=0;
ReadObj2.CurrentTime=0;
MCSval=1;
afterEcmT=[min(min(ecmTvar)), max(max(ecmTvar))];
while hasFrame(ReadObj)
    vidFrame = readFrame(ReadObj);
    aftervector{MCSval,1}=vidFrame;
    
    aftervector{MCSval,2}=rot90(reshape(ecmTvar(:,MCSval),grid,grid));
    
    vidFrame2 = readFrame(ReadObj2);
    aftervector{MCSval,3}=vidFrame2;
    MCSval=MCSval+1;
end


fullvector=cell(beforetime+aftertime,3);
for i=1:beforetime
    fullvector{i,1}=beforevector{i,1};
    fullvector{i,2}=beforevector{i,2};
    fullvector{i,3}=beforevector{i,3};
end
for i=beforetime+1:beforetime+aftertime
    fullvector{i,1}=aftervector{i-beforetime,1};
    fullvector{i,2}=aftervector{i-beforetime,2};
    fullvector{i,3}=aftervector{i-beforetime,3};
end
EcmtRange=[beforeEcmT; afterEcmT];
EcmtRange=[min(min(EcmtRange)), max(max(EcmtRange))];
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Track Cell phenotype
fullcellprops=cell(beforetime+aftertime,2);
for MCS=1:beforetime
    fullcellprops{MCS,1}=beforecelltype{MCS};
    fullcellprops{MCS,2}=beforecellpheno{MCS};
end

for MCS=1:aftertime
    fullcellprops{beforetime+MCS,1}=aftercelltype{MCS};
    fullcellprops{beforetime+MCS,2}=aftercellpheno{MCS};
end

Kpheno=zeros(beforetime+aftertime,3);
Spheno=zeros(beforetime+aftertime,3);
Tpheno=zeros(beforetime+aftertime,3);
for MCS=1:beforetime+aftertime
    for c=1:size(fullcellprops{MCS,1},1)
        if(fullcellprops{MCS,1}(c)==1) %Kupffer
            if(0 <= fullcellprops{MCS,2}(c) && fullcellprops{MCS,2}(c)< 0.48)
                Kpheno(MCS,1)=Kpheno(MCS,1)+1;
            elseif(0.48 <= fullcellprops{MCS,2}(c) && fullcellprops{MCS,2}(c)< 0.76)
                Kpheno(MCS,2)=Kpheno(MCS,2)+1;
            else %if(2/3 <= cellA && cellA <=1)
                Kpheno(MCS,3)=Kpheno(MCS,3)+1;
            end
            
        elseif(fullcellprops{MCS,1}(c)==2) %Stellate
            if(0 <= fullcellprops{MCS,2}(c) && fullcellprops{MCS,2}(c)< 0.48)
                Spheno(MCS,1)=Spheno(MCS,1)+1;
            elseif(0.48 <= fullcellprops{MCS,2}(c) && fullcellprops{MCS,2}(c)< 0.76)
                Spheno(MCS,2)=Spheno(MCS,2)+1;
            else %if(2/3 <= cellA && cellA <=1)
                Spheno(MCS,3)=Spheno(MCS,3)+1;
            end
        elseif(fullcellprops{MCS,1}(c)==3) %Tumor
            if(0 <= fullcellprops{MCS,2}(c) && fullcellprops{MCS,2}(c)< 0.48)
                Tpheno(MCS,1)=Tpheno(MCS,1)+1;
            elseif(0.48 <= fullcellprops{MCS,2}(c) && fullcellprops{MCS,2}(c)< 0.76)
                Tpheno(MCS,2)=Tpheno(MCS,2)+1;
            else %if(2/3 <= cellA && cellA <=1)
                Tpheno(MCS,3)=Tpheno(MCS,3)+1;
            end
        end
        
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numFrames = ceil(ReadObj.FrameRate*ReadObj.Duration);
f=figure;
colormap(jet);
ReadObj.CurrentTime=0;
ReadObj2.CurrentTime=0;
MCSval=1;
vidfile = VideoWriter(['BeforeandAfterTumorwithcells&EcmT&typeTest1_' before(end-6:end) '_' num2str(videocount)  '.mp4'],'MPEG-4');
open(vidfile);
for i=1:(beforetime+aftertime)
    
    subplot(1,3,1)
    imshow(fullvector{i,3}) % Cell Type
    
    subplot(1,3,2)
    imshow(fullvector{i,1}) % Cell Activation/Phenotype
    
    subplot(1,3,3) % ECM Concentration
    imagesc(fullvector{i,2});
    set(gca,'YColor','white');
    set(gca,'XColor','white');
    title(['EcmT at MCS: ' num2str(i)]);
    
    caxis(EcmtRange);
    colorbar
    
    
    
%     subplot(2,3,4) % Kupffer phenotype
%     plot(1:i,Kpheno(1:i,1)','LineWidth',2,'Color','#0072BD')
%     hold on
%     plot(1:i,Kpheno(1:i,2)','LineWidth',2,'Color','#EDB120')
%     hold on
%     plot(1:i,Kpheno(1:i,3)','LineWidth',2,'Color','#A2142F')
%     leg=legend('Less', 'Intermediate', 'Most');
%     hold off
%     title(leg,'Transformed');
%     xlabel('MCS')
%     ylabel('Number of cells')
%     ylim([0,max(max(Kpheno))])
%     xlim([1,beforetime+aftertime])
%     title('Kupffer Cells')
%     
%     subplot(2,3,5)% Stellate phenotype
%     plot(1:i,Spheno(1:i,1)','LineWidth',2,'Color','#0072BD')
%     hold on
%     plot(1:i,Spheno(1:i,2)','LineWidth',2,'Color','#EDB120')
%     hold on
%     plot(1:i,Spheno(1:i,3)','LineWidth',2,'Color','#A2142F')
%     leg=legend('Less', 'Intermediate', 'Most');
%     hold off
%     title(leg,'Transformed');
%     xlabel('MCS')
%     ylabel('Number of cells')
%     ylim([0,max(max(Spheno))])
%     xlim([1,beforetime+aftertime])
%     title(['Stellate Cells'])
%     
%     subplot(2,3,6)% Tumor phenotype
%     plot(1:i,Tpheno(1:i,1)','LineWidth',2,'Color','#0072BD')
%     hold on
%     plot(1:i,Tpheno(1:i,2)','LineWidth',2,'Color','#EDB120')
%     hold on
%     plot(1:i,Tpheno(1:i,3)','LineWidth',2,'Color','#A2142F')
%     leg=legend('Less', 'Intermediate', 'Most');
%     hold off
%     title(leg,'Transformed');
%     xlabel('MCS')
%     ylabel('Number of cells')
%     ylim([0,max(max(Tpheno))])
%     xlim([1,beforetime+aftertime])
%     title('Tumor Cells')
    
    %MCSval=MCSval+1;
    
    f.Position = [2,143,1912,475];
    F = getframe(f);
    [RGB] = frame2im(F);
    writeVideo(vidfile, RGB);
    pause(1/ReadObj.FrameRate);
    cla reset
end
close(vidfile)
close all
pause(10)

toc

%% Indiviual cell type
ReadObj = VideoReader(['KType_test_100p_501mcs_5tgfb__BeforeTumor_divf_0.5_param_14_scale_0.5ECMscaled0_1000006' '.mp4']);
ReadObj2 = VideoReader(['SType_test_100p_501mcs_5tgfb__BeforeTumor_divf_0.5_param_14_scale_0.5ECMscaled0_1000006'  '.mp4']);
ReadObj3 = VideoReader(['TType_test_100p_501mcs_5tgfb__BeforeTumor_divf_0.5_param_14_scale_0.5ECMscaled0_1000006' '.mp4']);

numFrames = ceil(ReadObj.FrameRate*ReadObj.Duration);
ReadObj.CurrentTime=0;
ReadObj2.CurrentTime=0;
ReadObj3.CurrentTime=0;
MCSval=1;
KSTvideo=cell(beforetime+aftertime,3);
while hasFrame(ReadObj)
    vidFrame = readFrame(ReadObj);
    subplot(1,3,1)
    KSTvideo{MCSval,1}=vidFrame;
    
    vidFrame2 = readFrame(ReadObj2);
    subplot(1,3,2)
     KSTvideo{MCSval,2}=vidFrame2;
     
    vidFrame3 = readFrame(ReadObj3);
    subplot(1,3,3)
     KSTvideo{MCSval,3}=vidFrame3;
    MCSval=MCSval+1;
    
    
end


ReadObj = VideoReader(['KType_test_100p_508mcs_0tgfb__AfterTumor_divf_0.75_param_14_scale_0.5ECMscaled0_1000006' '.mp4']);
ReadObj2 = VideoReader(['SType_test_100p_508mcs_0tgfb__AfterTumor_divf_0.75_param_14_scale_0.5ECMscaled0_1000006'  '.mp4']);
ReadObj3 = VideoReader(['TType_test_100p_508mcs_0tgfb__AfterTumor_divf_0.75_param_14_scale_0.5ECMscaled0_1000006' '.mp4']);

numFrames = ceil(ReadObj.FrameRate*ReadObj.Duration);
ReadObj.CurrentTime=0;
ReadObj2.CurrentTime=0;
ReadObj3.CurrentTime=0;
while hasFrame(ReadObj)
    vidFrame = readFrame(ReadObj);
    subplot(1,3,1)
    KSTvideo{MCSval,1}=vidFrame;
    
    vidFrame2 = readFrame(ReadObj2);
    subplot(1,3,2)
     KSTvideo{MCSval,2}=vidFrame2;
     
    vidFrame3 = readFrame(ReadObj3);
    subplot(1,3,3)
     KSTvideo{MCSval,3}=vidFrame3;
    MCSval=MCSval+1;
    
    
end

%%
f=figure;
vidfile = VideoWriter(['Celltype_Phenotype' before(end-6:end) '.mp4'],'MPEG-4');
open(vidfile);

for i=1:beforetime+aftertime
    subplot(1,3,1)
    imshow(KSTvideo{i,1})
    subplot(1,3,2)
    imshow(KSTvideo{i,2})
    subplot(1,3,3)
    imshow(KSTvideo{i,3})
f.Position = [2,453,1885,525];
    F = getframe(f);
    [RGB] = frame2im(F);
    writeVideo(vidfile, RGB);
    pause(1/ReadObj.FrameRate);

end
close(vidfile)
close all