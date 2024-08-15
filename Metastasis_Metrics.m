%% Repeated trials metrics
% 1. Before- 
% a- average kupffer and stellate phenotype and average ECM
% NOT NECESSARY- 2. After- total distance for first few seeded tumor cells and average tumor
%phenotype
% 2 After-
% a-ECM under tumor- and average tumor phenotype 
% 3.After-
% a- On adjacency Matrix- find tumor cells and total non-tumor cell
% connections for each cell over time
% b- find initially seeded tumor cells and calculate average of tumo cells
% from each other over time Growth Factor (?M)

% 1
X1Label='K more sensitive';
variables=11;
xv=[1:1:variables];
legendlabel=['0','0.5','1','1.5','2','2.5','3','3.5','4','4.5','5'];
vartick={'0','0.5','1','1.5','2','2.5','3','3.5','4','4.5','5'};
repeat=10;
%% Load IDVal
%%% For Kup/Stel Ratio 
% IDrange1=reshape(1000373:1:1000414,6,7);
% IDrange2=reshape(1000520:1:1000547,4,7);
% IDrange3=reshape(1000644:1:1000663,10,2);
% IDrange4=[IDrange1;IDrange2];
% IDrange=zeros(10,9);
% IDrange(:,1)=IDrange4(:,1); %K1
% IDrange(:,2)=IDrange3(:,1); %K0.9
% IDrange(:,3)=IDrange4(:,2); %K0.8
% IDrange(:,4)=IDrange3(:,2); %K0.7
% IDrange(:,5:9)=IDrange4(:,3:7);


%IDval=1000223;
IDrange=[1000664:1:1000773];
IDrange=reshape(IDrange,10,variables);
% IDrange=reshape(IDrange,6,variables);
% %IDrange=reshape(IDrange,10,variables);
%  IDval2=1000548:1000587;
% IDval2=reshape(IDval2,4,variables);
%  IDrange=[IDrange;IDval2];
% repeat=10;
% IDrange=[1000774:1:1000883];
% IDrange=reshape(IDrange,repeat,variables);
%variables=15;
%% Phase 1
Beforefilename=cell(repeat,variables);



%filename='test_100p_800mcs_5tgfb__BeforeTumor_divf_0.5_param_14_scale_0.5ECMscaled0_';

for i=1:size(IDrange,1)*size(IDrange,2)
    IDname=dir(['*BeforeTumor*' num2str(IDrange(i)) '.mat']);	
    Beforefilename{i}=[IDname.name];
end

BeforeStellate=zeros(repeat,variables);
BeforeKupffer=zeros(repeat,variables);
BeforeECM=zeros(repeat,variables);
BeforeCellNum=zeros(repeat,variables);

for varID=1:variables
for repeatID=1:repeat
tic
before=[Beforefilename{repeatID,varID}];
load([before]);

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
BeforeKupffer(repeatID,varID)=mean(phenotype{end}(pop_index{end}==1));
BeforeStellate(repeatID,varID)=mean(phenotype{end}(pop_index{end}==2));
BeforeCellNum(repeatID,varID)=size(csize{end},1);
BeforeECM(repeatID,varID)=mean(ecmTvar(:,end));

toc
end
end
%% phase 1- Visualize various cell phenotypes and ECM over time
f=figure
Y1Label='Mean ECM Concentration (?M)';
Y2Label='Mean Cell Activation';

yyaxis left
plot(xv,mean(BeforeECM),'LineWidth',2);
ylabel(Y1Label,'FontSize',20)
ylim([0 0.07])
yyaxis right
plot(xv,mean(BeforeKupffer),'LineWidth',2);
hold on
plot(xv,mean(BeforeStellate),'LineWidth',2);
ylabel(Y2Label,'FontSize',20);
xlim([xv(1),xv(end)])
ylim([0 1])
xlabel(X1Label,'FontSize',20);
legend('ECM','Kupffer','Stellate');
ax=gca;
ax.FontSize=16;
xticks(xv);
xticklabels(vartick);
%xticklabels({KSfile{1,2},KSfile{2,2},KSfile{3,2},KSfile{4,2},KSfile{5,2}});
%title('Phase 1- Before Tumor','FontSize',20);
f.Position = [799 49 1000 400];
beep
% NOT NEEDED 2- Tumor Migration Calculation

% Afterfilename=cell(repeat,variables);
% % filename='test_100p_1000mcs_0tgfb__AfterTumor_divf_0.75_param_14_scale_0.5ECMscaled0_';
% 
% 
% for i=1:size(IDrange,1)*size(IDrange,2)
%     IDname=dir(['*AfterTumor*' num2str(IDrange(i)) '.mat']);	
%     Afterfilename{i}=[IDname.name];
% end
% 
% AfterTumorPheno=zeros(repeat,variables);
% AfterTumorDist=cell(repeat,variables);
% AfterTumorCellType=cell(repeat,variables);
% 
% for varID=1:variables
% for repeatID=1:repeat
% tic
% after=[Afterfilename{repeatID,varID}];
% load(after)
% AfterTumorPheno(repeatID,varID)=mean(phenotype{end}(pop_index{1}==3));
% TumorID=find(pop_index{1}==3);
% cells=size(csize{1},1);
% XY=cell(cells,size(csize,2));
% DistStoE=zeros(cells,1);
% for i=1:cells
%     for j=1:size(csize,2)
%     loc=reshape(ctag(:,j),100,100);
%     [f1,f2]=find(loc == i);
%     coordinates=[round(mean(f2),0), round(mean(f1),0)];
%     XY{i,j}=coordinates;
%     end
% end
% 
% for i=1:cells 
%     for j=2:size(csize,2)
%         First=XY{i,j-1};
%         Second=XY{i,j};
%     DistStoE(i)=DistStoE(i) + sqrt((Second(1)-First(1))^2+(Second(2)-First(2))^2);
%     end
% end
% 
% 
% AfterTumorDist{repeatID,varID}=DistStoE;
% AfterTumorCellType{repeatID,varID}=TumorID;
% disp([num2str(repeatID) ',' num2str(varID)]);
% end
% end
% %%
% TumorCellDist=cell(1,variables);
% TumorBoxplot=double.empty;
% xBox=double.empty;
% for i=1:variables
%     for j=1:repeat
%         TumorCellDist{i}=[TumorCellDist{i} ; AfterTumorDist{j,i}(AfterTumorCellType{j,i})];
%     end
%     TumorBoxplot=[TumorBoxplot; TumorCellDist{i}];
%     xBox=[xBox; i*ones(i*repeat,1)];
%     %TumorCellDist{i}=mean(TumorCellDist{i});
% end
% Y1Label='Total Distance Travelled (pixels)';
% Y2Label='Average Cell Phenotype';
% X1Label='Number of Tumor Cells';
% figure
% yyaxis left
% boxplot(TumorBoxplot,xBox);
% ylabel(Y1Label);
% yyaxis right
% plot(mean(AfterTumorPheno));
% ylabel(Y2Label);
% xlabel(X1Label);


%% Phase 2 Find ECM of initial Tumor cells in After Tumor Run and tumor phenotype
% variables=11;
% repeat=6;
% 
% %IDrange=reshape(IDrange,repeat,variables);
% 
% IDval=1000223;
% IDrange=[IDval:1:IDval+65];
% %IDrange=reshape(IDrange,6,11);
% IDrange=reshape(IDrange,repeat,variables);

Afterfilename=cell(repeat,variables);
filename='test_100p_1000mcs_0tgfb__AfterTumor_divf_0.75_param_14_scale_0.5ECMscaled0_';

for i=1:size(IDrange,1)*size(IDrange,2)
IDname=dir(['*AfterTumor*' num2str(IDrange(i)) '.mat']);
Afterfilename{i}=[IDname.name];
end

AfterTumorPhenoMCS=cell(repeat,variables);
AfterECM=cell(repeat,variables);


RepeatedNonTumorConn=cell(repeat,variables);
RepeatedDistOverTime=cell(repeat,variables);
RepeatedTumorGrowth=cell(repeat,variables);

for varID=1:variables
for repeatID=1:repeat
tic
after=[Afterfilename{repeatID,varID}];
[repeatID,varID]
load([after]);

grid=100;
fulltime=size(csize,2);
beforetime=fulltime;
beforevector=cell(fulltime,3);
solTvar=zeros(grid*grid,fulltime);
solEvar=zeros(grid*grid,fulltime);
asmEvar=zeros(grid*grid,fulltime);
ecmTvar=zeros(grid*grid,fulltime);
TumorECM=zeros(1,fulltime);
TumorPhenoMCS=zeros(1,fulltime);

nontumorconnections=zeros(fulltime,1);
tumorgrowth=zeros(fulltime,1);
AverageDistOverTime=zeros(fulltime,1);

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
    
    TumorID=find(pop_index{1}==3);
    TumorECM(MCS)=mean(ecmTvar(ismember(maskG,TumorID)));
    TumorPhenoMCS(MCS)=mean(phenotype{MCS}(pop_index{1}==3));
    tumorgrowth(MCS)=numel(find(pop_index{MCS}==3));
    %Metastasis Metrics
       %1***************************************************
    Tumorcells=find(pop_index{MCS}==3); % find index of tumor cells
    adjacency1=graph(conn_mat{MCS}); % convert adjacency matrix to graph format
    Edges1=adjacency1.Edges(:,1); % only focus on edge pairs and not weights
    left1=Edges1.EndNodes(1:size(Edges1,1)); % left part of pair
    right1=Edges1.EndNodes(size(Edges1,1)+1:end); % right part of pair
    Tumorleft=(left1==Tumorcells); % if left vector has tumorcell index
    Tumorright=(right1==Tumorcells);% if right vector has tumorcell index
    pairs=sum(Tumorleft+Tumorright); % see if left and right both have tumor cells
    nontumor=sum((pairs==1)); % if only 1 of the pair is tumor then 1, else 0, so excludes tumor-tumor pairings
    nontumorconnections(MCS)=nontumor; 
    
    %2 ***************************************************
    loc1=reshape(ctag(:,MCS),grid,grid);
loc=rot90(loc1); % Rotates matrix to match video orientation

cells=find(pop_index{1}==3);
    cells2=size(cells,1);
cellxy=cell(cells2,1);
TumorDist=zeros(cells2,cells2);
for i=1:cells2
[f1,f2]=find(loc == cells(i));
    coordinates=[round(mean(f2),0), round(mean(f1),0)];
    cellxy{i}=coordinates;

end

for i=1:cells2
    for j=1:cells2
        First=cellxy{i};
        Second=cellxy{j};
    TumorDist(i,j)=sqrt((Second(1)-First(1))^2+(Second(2)-First(2))^2);

    end
end

DiagDistance=triu(TumorDist,1);
Avgdist=mean(nonzeros(DiagDistance));
AverageDistOverTime(MCS)=Avgdist;
   % ***************************************************
end
RepeatedNonTumorConn{repeatID,varID}=nontumorconnections;
RepeatedDistOverTime{repeatID,varID}=AverageDistOverTime;
RepeatedTumorGrowth{repeatID,varID}=tumorgrowth;
AfterECM{repeatID,varID}=TumorECM;
AfterTumorPhenoMCS{repeatID,varID}=TumorPhenoMCS;
toc
end

end


AveragedConnections=zeros(fulltime,variables);
AveragedTumorDistance=zeros(fulltime,variables);
for i=1:variables
    Avgval1=zeros(1000,1);
    Avgval2=zeros(1000,1);
    for j=1:repeat
        Avgval1(:)=Avgval1(:) + RepeatedNonTumorConn{j,i};
        Avgval2(:)=Avgval2(:) + RepeatedDistOverTime{j,i};
        
    end
    AveragedConnections(:,i)=Avgval1./repeat;
    AveragedTumorDistance(:,i)=Avgval2./repeat;
end

% Average ECM under Tumor and Average Tumor phenotype
c1=colormap(jet);
c2=round(linspace(1,size(c1,1),variables),0);
AfterTumorECMAverage=cell(1,variables);
AfterTumorPhenoAverage=cell(1,variables);
for i=1:variables
    AfterTumorECMAverage{i}=zeros(1,1000);
    AfterTumorPhenoAverage{i}=zeros(1,1000);
    for j=1:repeat
        AfterTumorECMAverage{i}=AfterTumorECMAverage{i} + AfterECM{j,i};
        AfterTumorPhenoAverage{i}=AfterTumorPhenoAverage{i} + AfterTumorPhenoMCS{j,i};
        
    end
    AfterTumorECMAverage{i}=AfterTumorECMAverage{i}/repeat;
    AfterTumorPhenoAverage{i}=AfterTumorPhenoAverage{i}/repeat;
   % TumorBoxplot=[TumorBoxplot; TumorCellDist{i}];
    %xBox=[xBox; i*ones(i*5,1)];
    %TumorCellDist{i}=mean(TumorCellDist{i});
end

Y1Label='Average ECM for Tumors';
Y2Label='Average Tumor Phenotype';
X1Axis='MCS';
% figure
% subplot(1,2,1)
% for i=1:variables
% plot(AfterTumorECMAverage{i},'Color',c1(c2(i),:));
% hold on
% end
% hold off
% ylabel(Y1Label);
% xlabel(X1Axis);
% subplot(1,2,2)
% for i=1:variables
% plot(AfterTumorPhenoAverage{i},'Color',c1(c2(i),:));
% hold on
% end
% hold off
% ylabel(Y2Label);
% xlabel(X1Axis);
% legend(vartick);

% Colorgraphs for Average ECM and Pheno for tumors
TumorPhenoAverageColorGraph=zeros(fulltime,variables);
TumorECMAverageColorGraph=zeros(fulltime,variables);
for i=1:size(AfterTumorECMAverage,2)
    TumorPhenoAverageColorGraph(:,i)=AfterTumorPhenoAverage{i};
    TumorECMAverageColorGraph(:,i)=AfterTumorECMAverage{i};
    
end
%% graphing phase 2
Y1Label='Mean ECM Concentration (?M)';
Y2Label='Mean Tumor Cell Activation';
micron=5; % 5 to adjust for size scaling factor that was changed from 2 to 0.5
af1=figure
subplot(1,2,1)
imagesc(TumorPhenoAverageColorGraph')
xlabel('Time (MCS)','FontSize',20)
ylabel(X1Label,'FontSize',20)
yticks([1:size(AveragedConnections,2)]);
yticklabels(vartick);
caxis([0 1])
c3=colormap(jet);
c4=colorbar;
%c4.Label.String=Y2Label;
title(Y2Label);
set(c4,'FontSize',16);
ax=gca;
ax.FontSize=16;
subplot(1,2,2)
imagesc(TumorECMAverageColorGraph')
caxis([0 0.07])
xlabel('Time (MCS)','FontSize',20)
ylabel(X1Label,'FontSize',20)
yticks([1:size(AveragedConnections,2)]);
yticklabels(vartick);
af1.Position = [799 49 1300 500];
ax=gca;
ax.FontSize=16;
c5=colormap(jet);
c6=colorbar;
title(Y1Label);
%c6.Label.String=Y1Label;
set(c6,'FontSize',16);
ax=gca;
ax.FontSize=16;

% Non-Tumor Connections and Distance between Tumor Cells 
af2=figure

subplot(1,2,1)
imagesc(AveragedConnections')
xlabel('Time (MCS)','FontSize',20)
ylabel(X1Label,'FontSize',20)
yticks([1:size(AveragedConnections,2)]);
yticklabels(vartick);
caxis([0 120]);
c3=colormap(jet);
c4=colorbar;
%c4.Label.String='Non-Tumor Cell Connections';
title({'Mean Non-Tumor','Cell Connections'});
set(c4,'FontSize',16);
ax=gca;
ax.FontSize=16;
subplot(1,2,2)
imagesc(AveragedTumorDistance'*micron)
xlabel('Time (MCS)','FontSize',20)
ylabel(X1Label,'FontSize',20)
yticks([1:size(AveragedConnections,2)]);
yticklabels(vartick);
caxis([0 45*micron]);
c5=colormap(jet);
c6=colorbar;
%c6.Label.String={'Average Distance between','Initial Tumor Cells'};
title({'Mean Distance between','Initial Tumor Cells (?m)'});
set(c6,'FontSize',16);
ax=gca;
ax.FontSize=16;
af2.Position = [799 49 1300 500];
beep

%% Graphing all 10 trials
variables=11;
fb=figure
%before
boxXvals=zeros(repeat,variables);
for b=1:repeat
boxXvals(b,:)=[0:0.5:5];
end
boxXvals=reshape(boxXvals,repeat*variables,1);
BeforeEcmbox=reshape(BeforeECM,repeat*variables,1);
BeforeKbox=reshape(BeforeKupffer,repeat*variables,1);
BeforeSbox=reshape(BeforeStellate,repeat*variables,1);
subplot(1,3,1)
boxplot(BeforeEcmbox,boxXvals);
ylim([0 0.07])
xlabel(X1Label,'FontSize',20);
ylabel('Mean ECM Concentration','FontSize',20);
ax=gca;
ax.FontSize=16;
subplot(1,3,2)
boxplot(BeforeKbox,boxXvals);
ylim([0 1])
xlabel(X1Label,'FontSize',20);
ylabel('Mean Kupffer Cell Activation','FontSize',20);
ax=gca;
ax.FontSize=16;
subplot(1,3,3)
boxplot(BeforeSbox,boxXvals);
ylim([0 1])
xlabel(X1Label,'FontSize',20);
ylabel('Mean Stellate Cell Activation','FontSize',20);
ax=gca;
ax.FontSize=16;
fb.Position = [6 -156 1890 774];

%%

%after
% RepeatedDistOverTime
% RepeatedNonTumorConn
% AfterECM
% AfterTumorPhenoMCS
fa=figure
var2graph=RepeatedDistOverTime;

maxvals=zeros(repeat,variables);
for v1=1:2:variables
    for v2=1:1:repeat
        maxvals(v2,v1)=max(var2graph{v2,v1});
    end
end
varlabel={'Average Distance between','Initially Seeded Tumor Cells'};
for sub=1:6
    subplot(2,3,sub);
    for trial=1:repeat
        plot(var2graph{trial,sub*2-1},'LineWidth',2);
        if (sub==4 || sub==5 || sub==6)
        xlabel('Time (MCS)','FontSize',20)
        end
        if (sub==1 || sub==4)
        ylabel(varlabel,'FontSize',20)
        end
        title([vartick{sub*2-1} ' ?M'],'FontSize',20);
        ax=gca;
ax.FontSize=16;
ylim([0 max(max(maxvals))])
        hold on
    end
    hold off
end
fa.Position = [16 -276 1882 894];
% 34,174,1882,671
%%

c1=colormap(jet);
c2=round(linspace(1,size(c1,1),variables),0);
figure
subplot(1,2,1)
for i=1:variables
plot(AveragedConnections(1:end,i),'Color',c1(c2(i),:))
hold on
end
hold off
xlabel('MCS')
ylabel('Number of Non-Tumor Connections')


subplot(1,2,2)

for i=1:variables
plot(AveragedTumorDistance(1:end,i),'Color',c1(c2(i),:))
hold on
end
xlabel('MCS')
ylabel('Average distance between initial tumor cells')

legend('0','0.5','1','1.5','2','2.5','3','3.5','4','4.5','5');

%% Tumor Growth Metrics quantified using .mat files
TumorGrowthAverage=cell(1,variables);
for i=1:variables
    TumorGrowthAverage{i}=zeros(1,1000);
    for j=1:repeat
        TumorGrowthAverage{i}=TumorGrowthAverage{i} + RepeatedTumorGrowth{j,i}';
    end
    TumorGrowthAverage{i}=TumorGrowthAverage{i}'/repeat;
end

%% Visualizing Tumor Growth
%µ
title1={'Number of Tumor Cells'};
vartick={'1','2','3','4','5','6','7','8','9','10','11','12','13','14'};
ylabel1='Time of Tumor Entry (100s)';
xlabel1='Time (MCS)';
growthparam=TumorGrowthAverageTumorEntryTimeParam;
growthcolormap=zeros(size(growthparam,2),1000);
for i=1:size(growthparam,2)
    growthcolormap(i,:)=growthparam{i}';
end
if growthmax<max(max(growthcolormap))
    growthmax=max(max(growthcolormap));
end
growthmax;
fa=figure;
imagesc(growthcolormap)
c5=colormap(jet);
c6=colorbar;
set(c6,'FontSize',16);
title(title1);
xlabel(xlabel1,'FontSize',20);
ylabel(ylabel1,'FontSize',20);
caxis([0 round(growthmax/10,0)*10]);
yticklabels(vartick);
ax=gca;
ax.FontSize=16;
fa.Position = [781   407   478   413];