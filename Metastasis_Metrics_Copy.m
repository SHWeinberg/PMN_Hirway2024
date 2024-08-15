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
% from each other over time

% 1
X1Label='Number of Tumor Cells';
xv=[1:1:10];
legendlabel=['1','2','3','4','5','6','7','8','9','10'];
vartick={'1','2','3','4','5','6','7','8','9','10'};

variables=10;
repeat=6;
IDval=1000415;
IDrange=[IDval:1:IDval+59];
%IDrange=reshape(IDrange,6,11);
IDrange=reshape(IDrange,repeat,variables);
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
% 1- Visualize various cell phenotypes and ECM over time
figure
Y1Label='Average ECM Concentration';
Y2Label='Average Cell Phenotype';

yyaxis left
plot(xv,mean(BeforeECM),'LineWidth',2);
ylabel(Y1Label)
yyaxis right
plot(xv,mean(BeforeKupffer),'LineWidth',2);
hold on
plot(xv,mean(BeforeStellate),'LineWidth',2);
ylabel(Y2Label);
xlim([xv(1),xv(end)])
ylim([0 1])
xlabel(X1Label);
legend('ECM','Kupffer','Stellate');
xticks(xv);
xticklabels(vartick);
%xticklabels({KSfile{1,2},KSfile{2,2},KSfile{3,2},KSfile{4,2},KSfile{5,2}});
title('Final- Before Tumor');
%% Graphing each repeated trial for variables to check 
figure
Y1Label='Average ECM Concentration';
Y2Label='Average Cell Phenotype';


subplot(1,3,1)
boxplot(BeforeECM);
ylabel(Y1Label)
%xlim([xv(1),xv(end)])
xlabel(X1Label);
%legend('ECM','Kupffer','Stellate');
xticks(xv);
xticklabels(vartick);
title('Before-ECM');
subplot(1,3,2)
boxplot(BeforeKupffer);
ylabel(Y2Label);
%xlim([xv(1),xv(end)])
ylim([0 1])
xlabel(X1Label);
%legend('ECM','Kupffer','Stellate');
xticks(xv);
xticklabels(vartick);
title('Before-Kupffer Pheno');
subplot(1,3,3)
boxplot(BeforeStellate);
ylabel(Y2Label);
%xlim([xv(1),xv(end)])
ylim([0 1])
xlabel(X1Label);
%legend('ECM','Kupffer','Stellate');
xticks(xv);
xticklabels(vartick);

%xticklabels({KSfile{1,2},KSfile{2,2},KSfile{3,2},KSfile{4,2},KSfile{5,2}});
title('Before-Stellate Pheno');

%% NOT NEEDED 2- Tumor Migration Calculation

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


%% 2 Find ECM of initial Tumor cells in After Tumor Run and tumor phenotype
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
figure
subplot(1,2,1)
imagesc(TumorPhenoAverageColorGraph')
xlabel('Time (MCS)')
ylabel(X1Label)
yticks([1:size(AveragedConnections,2)]);
yticklabels(vartick);
c3=colormap(jet);
c4=colorbar;
c4.Label.String=Y2Label;
subplot(1,2,2)
imagesc(TumorECMAverageColorGraph')
xlabel('Time (MCS)')
ylabel(X1Label)
yticks([1:size(AveragedConnections,2)]);
yticklabels(vartick);
c5=colormap(jet);
c6=colorbar;
c6.Label.String=Y1Label;


% Non-Tumor Connections and Distance between Tumor Cells 
figure
subplot(1,2,1)
imagesc(AveragedConnections')
xlabel('Time (MCS)')
ylabel(X1Label)
yticks([1:size(AveragedConnections,2)]);
yticklabels(vartick);
c3=colormap(jet);
c4=colorbar;
c4.Label.String='Non-Tumor Connections';
subplot(1,2,2)
imagesc(AveragedTumorDistance')
xlabel('Time (MCS)')
ylabel(X1Label)
yticks([1:size(AveragedConnections,2)]);
yticklabels(vartick);
c5=colormap(jet);
c6=colorbar;
c6.Label.String='Average Distance between Initial Tumor Cells';
%% Graphing each repeated trials
xv2=[1 3 5 7 9 10];
figure
%subplot (2,3,1)
for j=1:size(xv2,2)
    subplot(2,3,j)
for i=1:repeat
    plot(AfterTumorPhenoMCS{i,xv2(j)},'LineWidth',2)
    hold on
end
ylabel('Tumor Phenotype');
title([vartick{xv2(j)} ' Tumor cells'])
xlabel('MCS');
%legend('1','2','3','4','5','6');
end
%**************************************
figure
subplot (2,3,1)
for j=1:size(xv2,2)
    subplot(2,3,j)
for i=1:repeat
    plot(AfterECM{i,xv2(j)},'LineWidth',2)
    hold on
end
ylabel('ECM for Tumors');
title([vartick{xv2(j)} ' Tumor cells'])
xlabel('MCS');
%legend('1','2','3','4','5','6');
end
%**************************************
figure
subplot (2,3,1)
for j=1:size(xv2,2)
    subplot(2,3,j)
for i=1:repeat
    plot(RepeatedNonTumorConn{i,xv2(j)},'LineWidth',2)
    hold on
end
ylabel('Non-Tumor Connections');
title([vartick{xv2(j)} ' Tumor cells'])
xlabel('MCS');
%legend('1','2','3','4','5','6');
end
%**************************************
figure
subplot (2,3,1)
for j=1:size(xv2,2)
    subplot(2,3,j)
for i=1:repeat
    plot(RepeatedDistOverTime{i,xv2(j)},'LineWidth',2)
    hold on
end
hold off
ylabel('Avg Distance- Tumor Cells');
title([vartick{xv2(j)} ' Tumor cells'])
xlabel('MCS');
%legend('1','2','3','4','5','6');
end


%% Metastasis metrics

% variables=14;
% repeat=6;
% IDval=1000289;
% IDrange=[IDval:1:IDval+83];
% %IDrange=reshape(IDrange,6,11);
% IDrange=reshape(IDrange,repeat,variables);
%1



Afterfilename=cell(repeat,variables);
%filename='test_100p_800mcs_5tgfb__BeforeTumor_divf_0.5_param_14_scale_0.5ECMscaled0_';

for i=1:size(IDrange,1)*size(IDrange,2)
    IDname=dir(['*AfterTumor*' num2str(IDrange(i)) '.mat']);	
    Afterfilename{i}=[IDname.name];
end
dim=100;
% for ID1=1:variables
%     for ID2=1:repeat
%     [ID2,ID1]
% load (Afterfilename{ID2,ID1});
% 
% 
% for MCS=1:fulltime
%  
% 
% 
% 
%     end
% end

%



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
legend(legendlabel);

%legend('100','200','300','400','500','600','700','800','900','1000','1100','1200','1300','1400');

