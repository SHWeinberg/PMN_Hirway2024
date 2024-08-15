% Averaging over several trials
% Get names of test.mat files in folder
% sort by comment keyword
% sort by tgfb value and parameter
%%
% For whatever parameter
tic
temp_cell=cell(1,11); %Create temporary cell

filetype='*PMN_scf0.5.mat';
files1=GetNames(filetype);
files1=files1';

files2=sort(files1);
Index0 = find(contains(files2,'0tgfb'));
Index5 = find(contains(files2,'5tgfb'));
files0=files2(Index0);
files5=files2(Index5);

parameter='tgfb'; %'kDAE.scale1.05';

files0p=files0(find(contains(files0,parameter)));
files5p=files5(find(contains(files5,parameter)));

temp_cell{1,11}=parameter; % parameter name
temp_cell{1,1}=0; % 0 tgfb
temp_cell{1,2}=files0p; % filenames for 0 tgfb
temp_cell{1,6}=5; % 5 tgfb
temp_cell{1,7}=files5p; % filenames for 5 tgfb
toc


%% Loading files

tic
dim=100;

sizeAvg0=zeros(1,1200);
phenoAvg0=zeros(1,1200);
jforceAvg0=zeros(1,1200);
confAvg0=zeros(1,1200);

sizeAvg5=zeros(1,1200);
phenoAvg5=zeros(1,1200);
jforceAvg5=zeros(1,1200);
confAvg5=zeros(1,1200);
grid=100+100+98+98;

for sim=1:5
    
    load(files0p{sim}); %load 0 file (sim)
    for tp=1:1200
       % disp(['0 TGFB time: ' num2str(tp)]);
        
        cells=size(csize{1,tp},1);
        loc1=reshape(ctag(:,tp),dim,dim);
        loc=rot90(loc1); % Rotates matrix to match video orientation
        size(loc);
        finalsize=csize{1,tp};
        finalp=phenotype{1,tp};
        s1=zeros(dim,dim);
        p1=zeros(dim,dim);
        % s2=zeros(dim,dim);
        % p2=zeros(dim,dim);
        
        %Centroid mapped over dim x dim matrix
        %     for i=1:cells
        %         [f1,f2]=find(loc == i);
        %         xc=round(mean(f2),0);
        %         yc=round(mean(f1),0);
        %         s1(xc,yc)=finalsize(i);
        %         p1(xc,yc)=finalp(i);
        %         s2(f1,f2)=finalsize(i);
        %         p2(f1,f2)=finalp(i);
        %     end
        %     s1=s1';%transposing makes it match video - b1 matrix with cell size at centroid
        %     p1=p1'; %p1 matrix for phenotype
        
        % Overall Average
        
        s2=mean(finalsize);
        p2=mean(finalp);
        conf=sum(finalsize)/(dim*dim-grid);
        
        
        
        sizeAvg0(:,tp)=sizeAvg0(:,tp)+s2;
        phenoAvg0(:,tp)=phenoAvg0(:,tp)+p2;
        confAvg0(:,tp)=confAvg0(:,tp)+conf;
    end
    
    
    %5 tgfb
    load(files5p{sim}); % load 5tgfb file (sim)
    for tp=1:1200
       % disp(['5 TGFB time: ' num2str(tp)]);
        cells=size(csize{1,tp},1);
        loc1=reshape(ctag(:,tp),dim,dim);
        loc=rot90(loc1); % Rotates matrix to match video orientation
        size(loc);
        finalsize=csize{1,tp};
        finalp=phenotype{1,tp};
        s1=zeros(dim,dim);
        p1=zeros(dim,dim);
        % s2=zeros(dim,dim);
        %p2=zeros(dim,dim);
        
        
        %Centroid mapped over dim x dim matrix
        %     for i=1:cells
        %         [f1,f2]=find(loc == i);
        %         xc=round(mean(f2),0);
        %         yc=round(mean(f1),0);
        %         s1(xc,yc)=finalsize(i);
        %         p1(xc,yc)=finalp(i);
        %         s2(f1,f2)=finalsize(i);
        %         p2(f1,f2)=finalp(i);
        %     end
        %     s1=s1';%transposing makes it match video - b1 matrix with cell size at centroid
        %     p1=p1'; %p1 matrix for phenotype
        
        
        s2=mean(finalsize);
        p2=mean(finalp);
        conf=sum(finalsize)/(dim*dim-grid);
        
        sizeAvg5(:,tp)=sizeAvg5(:,tp)+s2;
        phenoAvg5(:,tp)=phenoAvg5(:,tp)+p2;
        confAvg5(:,tp)=confAvg5(:,tp)+conf;
    end
end
%Averaging
sizeAvg0=sizeAvg0/5;
phenoAvg0=phenoAvg0/5;
confAvg0=confAvg0/5;

sizeAvg5=sizeAvg5/5;
phenoAvg5=phenoAvg5/5;
confAvg5=confAvg5/5;


temp_cell{1,3}=phenoAvg0;
temp_cell{1,4}=sizeAvg0;
temp_cell{1,5}=confAvg0;

temp_cell{1,8}=phenoAvg5;
temp_cell{1,9}=sizeAvg5;
temp_cell{1,10}=confAvg5;



% Junctional Forces
% ctag_mat=reshape(ctag(:,tp),dim,dim);
%     [xjunc_cent, yjunc_cent, jx_cent, jy_cent, ...
%         xjunc_vox, yjunc_vox, jx_vox_norm, jy_vox_norm, jx_vox, jy_vox,...
%         jx_mat, jy_mat, neighbors_mat, Aconn] ...
%         = calc_junc_from_ctag(ctag_mat, dim, dim);
toc
beep
%% Store temp cell in the proper variable
for i=1:11
jhalf_var2{4,i}=temp_cell{1,i};
end

%'$e_1/e_\mathrm{in}$'
%% Graphing 
colorj=jet;
graph_var=tmax_var;
basename='Jhalf9574';
baseval='1';
legtitle='T_{max}';
delx=2.5*2.5;
fa=24;
fb=fa+2;
fc=fb+2;

colormin=min(cell2mat(graph_var(:,12)));
colormax=max(cell2mat(graph_var(:,12)));
colorrange=colormin:(colormax-colormin)/63:colormax;
[minValue,closestIndex1] = min(abs(colorrange-graph_var{1,12}));
[minValue,closestIndex2] = min(abs(colorrange-graph_var{2,12}));
[minValue,closestIndex3] = min(abs(colorrange-1));
[minValue,closestIndex4] = min(abs(colorrange-graph_var{3,12}));
[minValue,closestIndex5] = min(abs(colorrange-graph_var{4,12}));


%0 tgfb
figure
tp=1:1200;
subplot(2,3,1)
plot(graph_var{1,3},'LineWidth',2,'Color',colorj(closestIndex1,:));
hold on
plot(graph_var{2,3},'LineWidth',2,'Color',colorj(closestIndex2,:));
hold on
plot(base_var{1,3},'LineWidth',2,'Color',colorj(closestIndex3,:));
hold on
plot(graph_var{3,3},'LineWidth',2,'Color',colorj(closestIndex4,:));
hold on
plot(graph_var{4,3},'LineWidth',2,'Color',colorj(closestIndex5,:));
ylim([0 1]);
a = get(gca,'XTick');
set(gca,'XTick',a,'fontsize',fa-2)
b = get(gca,'YTick');
set(gca,'YTick',b,'fontsize',fa-2)

%xlabel('Time (MCS)','FontSize',fa);
xticks([0 400 800 1200])
xticklabels({'0','400','800','1200'})
yticks([0 0.25 0.5 0.75 1])
yticklabels({'0','0.25','0.5','0.75','1'})
ylabel('Cell Phenotype','FontSize',fa);
set(gca,'box','off') 
title('Phenotype','FontSize',fb,'FontWeight','normal');
%legend(graph_var{1,11},graph_var{2,11},basename,graph_var{3,11},graph_var{4,11});

sizev1=[graph_var{1,4}*delx,graph_var{2,4}*delx,base_var{1,4}*delx,graph_var{3,4}*delx,graph_var{4,4}*delx];
maxv1=max(max(sizev1));
maxrange1=maxv1+(100-mod(maxv1,100));
subplot(2,3,2)
plot(graph_var{1,4}*delx,'LineWidth',2,'Color',colorj(closestIndex1,:));
hold on
plot(graph_var{2,4}*delx,'LineWidth',2,'Color',colorj(closestIndex2,:));
hold on
plot(base_var{1,4}*delx,'LineWidth',2,'Color',colorj(closestIndex3,:));
hold on
plot(graph_var{3,4}*delx,'LineWidth',2,'Color',colorj(closestIndex4,:));
hold on
plot(graph_var{4,4}*delx,'LineWidth',2,'Color',colorj(closestIndex5,:));
ylabel('Cell Area (?m^2)','FontSize',fa);
%xlabel('Time (MCS)','FontSize',fa);
xticks([0 400 800 1200])
xticklabels({'0','400','800','1200'})
yticks([0 maxrange1/4 maxrange1*2/4 maxrange1*3/4 maxrange1])
yticklabels({num2str(0) num2str(maxrange1/4) num2str(maxrange1*2/4) num2str(maxrange1*3/4) num2str(maxrange1)})
title('Cell Size','FontSize',fb,'FontWeight','normal');
a = get(gca,'XTick');
set(gca,'XTick',a,'fontsize',fa-2)
b = get(gca,'YTick');
set(gca,'YTick',b,'fontsize',fa-2)
set(gca,'box','off') 
leg=legend(num2str(graph_var{1,12}),num2str(graph_var{2,12}),baseval,num2str(graph_var{3,12}),num2str(graph_var{4,12}),'FontSize',fa-2 ,'Interpreter', 'tex');
title(leg,legtitle);

subplot(2,3,3)
plot(graph_var{1,5},'LineWidth',2,'Color',colorj(closestIndex1,:));
hold on
plot(graph_var{2,5},'LineWidth',2,'Color',colorj(closestIndex2,:));
hold on
plot(base_var{1,5},'LineWidth',2,'Color',colorj(closestIndex3,:));
hold on
plot(graph_var{3,5},'LineWidth',2,'Color',colorj(closestIndex4,:));
hold on
plot(graph_var{4,5},'LineWidth',2,'Color',colorj(closestIndex5,:));
%xlabel('Time (MCS)','FontSize',fa);
xticks([0 400 800 1200])
xticklabels({'0','400','800','1200'})
yticks([0 0.25 0.5 0.75 1])
yticklabels({'0','0.25','0.5','0.75','1'})
ylabel('Grid Confluence','FontSize',fa);
title('Confluence','FontSize',fb,'FontWeight','normal');
a = get(gca,'XTick');
set(gca,'XTick',a,'fontsize',fa-2)
b = get(gca,'YTick');
set(gca,'YTick',b,'fontsize',fa-2)
set(gca,'box','off') 
caxis([colormin colormax])
colormap(colorj)

%legend(graph_var{1,11},graph_var{2,11},basename,graph_var{3,11},graph_var{4,11});
%sgtitle('Measurements for 0 TGF-?','FontSize',fc);

% 5 Tgfb

tp=1:1200;
subplot(2,3,4)
plot(graph_var{1,8},'LineWidth',2,'Color',colorj(closestIndex1,:));
hold on
plot(graph_var{2,8},'LineWidth',2,'Color',colorj(closestIndex2,:));
hold on
plot(base_var{1,8},'LineWidth',2,'Color',colorj(closestIndex3,:));
hold on
plot(graph_var{3,8},'LineWidth',2,'Color',colorj(closestIndex4,:));
hold on
plot(graph_var{4,8},'LineWidth',2,'Color',colorj(closestIndex5,:));
ylim([0 1]);
a = get(gca,'XTick');
set(gca,'XTick',a,'fontsize',fa-2)
b = get(gca,'YTick');
set(gca,'YTick',b,'fontsize',fa-2)
xlabel('Time (MCS)','FontSize',fa);
xticks([0 400 800 1200])
xticklabels({'0','400','800','1200'})
yticks([0 0.25 0.5 0.75 1])
yticklabels({'0','0.25','0.5','0.75','1'})
ylabel('Cell Phenotype','FontSize',fa);
set(gca,'box','off') 
%title('Phenotype','FontSize',fb);
%legend(graph_var{1,11},graph_var{2,11},basename,graph_var{3,11},graph_var{4,11});

subplot(2,3,5)
sizev2=[graph_var{1,9}*delx,graph_var{2,9}*delx,base_var{1,9}*delx,graph_var{3,9}*delx,graph_var{4,9}*delx];
maxv2=max(max(sizev2));
maxrange2=maxv2+(100-mod(maxv2,100));
plot(graph_var{1,9}*delx,'LineWidth',2,'Color',colorj(closestIndex1,:));
hold on
plot(graph_var{2,9}*delx,'LineWidth',2,'Color',colorj(closestIndex2,:));
hold on
plot(base_var{1,9}*delx,'LineWidth',2,'Color',colorj(closestIndex3,:));
hold on
plot(graph_var{3,9}*delx,'LineWidth',2,'Color',colorj(closestIndex4,:));
hold on
plot(graph_var{4,9}*delx,'LineWidth',2,'Color',colorj(closestIndex5,:));
ylabel('Cell Area (?m^2)','FontSize',fa);
xlabel('Time (MCS)','FontSize',fa);
xticks([0 400 800 1200])
xticklabels({'0','400','800','1200'})
yticks([0 maxrange2/4 maxrange2*2/4 maxrange2*3/4 maxrange2])
yticklabels({num2str(0) num2str(maxrange2/4) num2str(maxrange2*2/4) num2str(maxrange2*3/4) num2str(maxrange2)})

a = get(gca,'XTick');
set(gca,'XTick',a,'fontsize',fa-2)
b = get(gca,'YTick');
set(gca,'YTick',b,'fontsize',fa-2)
set(gca,'box','off') 
title('5 TGF-?','FontSize',fc,'FontWeight','bold');
leg=legend(num2str(graph_var{1,12}),num2str(graph_var{2,12}),baseval,num2str(graph_var{3,12}),num2str(graph_var{4,12}),'FontSize',fa-2, 'Interpreter', 'tex');
title(leg,legtitle);

subplot(2,3,6)
plot(graph_var{1,10},'LineWidth',2,'Color',colorj(closestIndex1,:));
hold on
plot(graph_var{2,10},'LineWidth',2,'Color',colorj(closestIndex2,:));
hold on
plot(base_var{1,10},'LineWidth',2,'Color',colorj(closestIndex3,:));
hold on
plot(graph_var{3,10},'LineWidth',2,'Color',colorj(closestIndex4,:));
hold on
plot(graph_var{4,10},'LineWidth',2,'Color',colorj(closestIndex5,:));
xlabel('Time (MCS)','FontSize',fa);
xticks([0 400 800 1200])
xticklabels({'0','400','800','1200'})
yticks([0 0.25 0.5 0.75 1])
yticklabels({'0','0.25','0.5','0.75','1'})
ylabel('Grid Confluence','FontSize',fa);
a = get(gca,'XTick');
set(gca,'XTick',a,'fontsize',fa-2)
b = get(gca,'YTick');
set(gca,'YTick',b,'fontsize',fa-2)
set(gca,'box','off') 
%title('Confluence','FontSize',fb);

caxis([colormin colormax])
colormap(colorj)


%legend(graph_var{1,11},graph_var{2,11},basename,graph_var{3,11},graph_var{4,11});
sgtitle('0 TGF-?','FontSize',fc,'FontWeight','bold');
%%
% For TGFB or other individual analysis
tic
temp_cell=cell(1,6); %Create temporary cell

filetype='*TgfbRR.mat';
files1=GetNames(filetype);
files1=files1';

files2=sort(files1);
% Index0 = find(contains(files2,'0tgfb'));
% Index5 = find(contains(files2,'5tgfb'));
% files0=files2(Index0);
% files5=files2(Index5);

parameter='2.5tgfb'; %'kDAE.scale1.05';

filestgfb=files2(find(contains(files2,parameter)));
% files5p=files5(find(contains(files5,parameter)));
simsize=size(filestgfb,1);
temp_var=cell(simsize,6);
%tgfb_var{1,6}=parameter; % parameter name
%tgfb_var{1,1}=0; % 0 tgfb
%tgfb_var{1,2}=filestgfb; % filenames for 0 tgfb

toc


tic

dim=100;

sizeAvg0=zeros(1,1200);
phenoAvg0=zeros(1,1200);
jforceAvg0=zeros(1,1200);
confAvg0=zeros(1,1200);

grid=100+100+98+98;
%%
tic
   temp_var{1,2}=filestgfb([1,2,3,4,5]);
for sim=[1,2,3,4,5]
    
    load(filestgfb{sim}); %load 0 file (sim)
 
    for tp=1:1200
       % disp(['0 TGFB time: ' num2str(tp)]);
        
        cells=size(csize{1,tp},1);
        loc1=reshape(ctag(:,tp),dim,dim);
        loc=rot90(loc1); % Rotates matrix to match video orientation
        size(loc);
        finalsize=csize{1,tp};
        finalp=phenotype{1,tp};
        s1=zeros(dim,dim);
        p1=zeros(dim,dim);
        % s2=zeros(dim,dim);
        % p2=zeros(dim,dim);
        
        %Centroid mapped over dim x dim matrix
        %     for i=1:cells
        %         [f1,f2]=find(loc == i);
        %         xc=round(mean(f2),0);
        %         yc=round(mean(f1),0);
        %         s1(xc,yc)=finalsize(i);
        %         p1(xc,yc)=finalp(i);
        %         s2(f1,f2)=finalsize(i);
        %         p2(f1,f2)=finalp(i);
        %     end
        %     s1=s1';%transposing makes it match video - b1 matrix with cell size at centroid
        %     p1=p1'; %p1 matrix for phenotype
        
        % Overall Average
        
        s2=mean(finalsize);
        p2=mean(finalp);
        conf=sum(finalsize)/(dim*dim-grid);
        
        
        
        sizeAvg0(:,tp)=sizeAvg0(:,tp)+s2;
        phenoAvg0(:,tp)=phenoAvg0(:,tp)+p2;
        confAvg0(:,tp)=confAvg0(:,tp)+conf;
    end
    
end


%Averaging
sizeAvg0=sizeAvg0/5;
phenoAvg0=phenoAvg0/5;
confAvg0=confAvg0/5;

    temp_var{1,3}=phenoAvg0;
    temp_var{1,4}=sizeAvg0;
    temp_var{1,5}=confAvg0;

% Junctional Forces
% ctag_mat=reshape(ctag(:,tp),dim,dim);
%     [xjunc_cent, yjunc_cent, jx_cent, jy_cent, ...
%         xjunc_vox, yjunc_vox, jx_vox_norm, jy_vox_norm, jx_vox, jy_vox,...
%         jx_mat, jy_mat, neighbors_mat, Aconn] ...
%         = calc_junc_from_ctag(ctag_mat, dim, dim);
toc
beep
%%  adding original 0 tgfb to RR variable
rowi=4;
ka_varRR{rowi,1}=ka_var{rowi,1};
ka_varRR{rowi,2}=ka_var{rowi,2};
ka_varRR{rowi,3}=ka_var{rowi,3};
ka_varRR{rowi,4}=ka_var{rowi,4};
ka_varRR{rowi,5}=ka_var{rowi,5};
ka_varRR{rowi,12}=ka_var{rowi,12};
%% Storing temp variable

ka_varRR{rowi,6}=5;
ka_varRR{rowi,7}=temp_var{1,2};
ka_varRR{rowi,8}=temp_var{1,3};
ka_varRR{rowi,9}=temp_var{1,4};
ka_varRR{rowi,10}=temp_var{1,5};

%%
tgfb_varRR{4,1}=base_varRR{1,6};
tgfb_varRR{4,2}=base_varRR{1,7};
tgfb_varRR{4,3}=base_varRR{1,8};
tgfb_varRR{4,4}=base_varRR{1,9};
tgfb_varRR{4,5}=base_varRR{1,10};

%% comparing 5 and 1.5 tgfb
figure
subplot(1,3,1)
plot(tgfb_var2{1,3})
hold on
plot(tgfb_var2{2,3})
hold on
plot(temp_var{1,3})
hold on
plot(tgfb_var2{3,3})
hold on
plot(tgfb_var2{4,3})
xlabel('Time (MCS)')
ylabel('Phenotype')
subplot(1,3,2)
plot(tgfb_var2{1,4}*delx)
hold on
plot(tgfb_var2{2,4}*delx)
hold on
plot(temp_var{1,4}*delx)
hold on
plot(tgfb_var2{3,4}*delx)
hold on
plot(tgfb_var2{4,4}*delx)
xlabel('Time (MCS)')
ylabel('Cell Size')
legend('0','1','1.8','2.5','5');
subplot(1,3,3)
plot(tgfb_var2{1,5})
hold on
plot(tgfb_var2{2,5})
hold on
plot(temp_var{1,5})
hold on
plot(tgfb_var2{3,5})
hold on
plot(tgfb_var2{4,5})
xlabel('Time (MCS)')
ylabel('Confluence')









%% *********** Visualize 3x3 Graph of pheno, cell size, confluence and steady states for 0 and 5 tgfb ***************************
figure
colorj=turbo;
graph_var=kDAE_varRR;
basename='Jhalf9574';
baseval='1';
legtitle='kd_{AE}';
delx=2.5*2.5;
fa=22;
fb=fa+2;
fc=fb+2;

colormin=min(cell2mat(graph_var(:,12)));
colormax=max(cell2mat(graph_var(:,12)));
colorrange=colormin:(colormax-colormin)/255:colormax;
[minValue,closestIndex1] = min(abs(colorrange-graph_var{1,12}));
[minValue,closestIndex2] = min(abs(colorrange-graph_var{2,12}));
[minValue,closestIndex3] = min(abs(colorrange-1));
[minValue,closestIndex4] = min(abs(colorrange-graph_var{3,12}));
[minValue,closestIndex5] = min(abs(colorrange-graph_var{4,12}));


%1- 0 tgfb control pheno
tp=1:1200;
subplot(3,3,1)
plot(graph_var{1,3},'LineWidth',2,'Color',colorj(closestIndex1,:));
hold on
plot(graph_var{2,3},'LineWidth',2,'Color',colorj(closestIndex2,:));
hold on
plot(base_varRR{1,3},'LineWidth',2,'Color',colorj(closestIndex3,:));
hold on
plot(graph_var{3,3},'LineWidth',2,'Color',colorj(closestIndex4,:));
hold on
plot(graph_var{4,3},'LineWidth',2,'Color',colorj(closestIndex5,:));
ylim([0 1]);
a = get(gca,'XTick');
set(gca,'XTick',a,'fontsize',fa-2)
b = get(gca,'YTick');
set(gca,'YTick',b,'fontsize',fa-2)

%xlabel('Time (MCS)','FontSize',fa);
xticks([0 400 800 1200])
xticklabels({'0','400','800','1200'})
yticks([0 0.25 0.5 0.75 1])
yticklabels({'0','0.25','0.5','0.75','1'})
ylabel({'Cell','Phenotype'},'FontSize',fa);
set(gca,'box','off') 
title('Phenotype','FontSize',fa,'FontWeight','normal');
%legend(graph_var{1,11},graph_var{2,11},basename,graph_var{3,11},graph_var{4,11});

% 2- 0 tgfb control Cell Size
sizev1=[graph_var{1,4}*delx,graph_var{2,4}*delx,base_varRR{1,4}*delx,graph_var{3,4}*delx,graph_var{4,4}*delx];
sizev2=[graph_var{1,9}*delx,graph_var{2,9}*delx,base_varRR{1,9}*delx,graph_var{3,9}*delx,graph_var{4,9}*delx];
maxv1=max(max([sizev1 sizev2]));
maxrange1=maxv1+(100-mod(maxv1,100));

subplot(3,3,2)
plot(graph_var{1,4}*delx,'LineWidth',2,'Color',colorj(closestIndex1,:));
hold on
plot(graph_var{2,4}*delx,'LineWidth',2,'Color',colorj(closestIndex2,:));
hold on
plot(base_varRR{1,4}*delx,'LineWidth',2,'Color',colorj(closestIndex3,:));
hold on
plot(graph_var{3,4}*delx,'LineWidth',2,'Color',colorj(closestIndex4,:));
hold on
plot(graph_var{4,4}*delx,'LineWidth',2,'Color',colorj(closestIndex5,:));
ylabel({'Cell Area (μm^2)'},'FontSize',fa);
%xlabel('Time (MCS)','FontSize',fa);
xticks([0 400 800 1200])
xticklabels({'0','400','800','1200'})
yticks([300 300+(maxrange1-300)/4 300+2*(maxrange1-300)/4  300+3*(maxrange1-300)/4 maxrange1])
yticklabels({num2str(300) num2str(300+(maxrange1-300)/4) num2str(300+2*(maxrange1-300)/4) num2str(300+3*(maxrange1-300)/4) num2str(maxrange1)})
ylim([300 maxrange1])
title('Cell Size','FontSize',fa,'FontWeight','normal');
a = get(gca,'XTick');
set(gca,'XTick',a,'fontsize',fa-2)
b = get(gca,'YTick');
set(gca,'YTick',b,'fontsize',fa-2)
set(gca,'box','off') 
leg=legend(num2str(graph_var{1,12}),num2str(graph_var{2,12}),baseval,num2str(graph_var{3,12}),num2str(graph_var{4,12}),'FontSize',fa-6 ,'Interpreter', 'tex');
title(leg,legtitle);

% 3- 0 tgfb confluence
subplot(3,3,3)
plot(graph_var{1,5},'LineWidth',2,'Color',colorj(closestIndex1,:));
hold on
plot(graph_var{2,5},'LineWidth',2,'Color',colorj(closestIndex2,:));
hold on
plot(base_varRR{1,5},'LineWidth',2,'Color',colorj(closestIndex3,:));
hold on
plot(graph_var{3,5},'LineWidth',2,'Color',colorj(closestIndex4,:));
hold on
plot(graph_var{4,5},'LineWidth',2,'Color',colorj(closestIndex5,:));
%xlabel('Time (MCS)','FontSize',fa);
xticks([0 400 800 1200])
xticklabels({'0','400','800','1200'})
yticks([0 0.25 0.5 0.75 1])
yticklabels({'0','0.25','0.5','0.75','1'})
ylabel({'Grid','Confluence'},'FontSize',fa);
title('Confluence','FontSize',fa,'FontWeight','normal');
a = get(gca,'XTick');
set(gca,'XTick',a,'fontsize',fa-2)
b = get(gca,'YTick');
set(gca,'YTick',b,'fontsize',fa-2)
set(gca,'box','off') 
caxis([colormin colormax])
colormap(colorj)

%legend(graph_var{1,11},graph_var{2,11},basename,graph_var{3,11},graph_var{4,11});
%sgtitle('Measurements for 0 TGF-?','FontSize',fc);

% 4- 5 Tgfb treated pheno *****************************

tp=1:1200;
subplot(3,3,4)
plot(graph_var{1,8},'LineWidth',2,'Color',colorj(closestIndex1,:));
hold on
plot(graph_var{2,8},'LineWidth',2,'Color',colorj(closestIndex2,:));
hold on
plot(base_varRR{1,8},'LineWidth',2,'Color',colorj(closestIndex3,:));
hold on
plot(graph_var{3,8},'LineWidth',2,'Color',colorj(closestIndex4,:));
hold on
plot(graph_var{4,8},'LineWidth',2,'Color',colorj(closestIndex5,:));
ylim([0 1]);
a = get(gca,'XTick');
set(gca,'XTick',a,'fontsize',fa-2)
b = get(gca,'YTick');
set(gca,'YTick',b,'fontsize',fa-2)
xlabel('Time (MCS)','FontSize',fa);
xticks([0 400 800 1200])
xticklabels({'0','400','800','1200'})
yticks([0 0.25 0.5 0.75 1])
yticklabels({'0','0.25','0.5','0.75','1'})
ylabel({'Cell','Phenotype'},'FontSize',fa);
set(gca,'box','off') 
%title('Phenotype','FontSize',fb);
%legend(graph_var{1,11},graph_var{2,11},basename,graph_var{3,11},graph_var{4,11});

% 5- 5 TGFB Treated Cell Size *****************************
subplot(3,3,5)
sizev2=[graph_var{1,9}*delx,graph_var{2,9}*delx,base_varRR{1,9}*delx,graph_var{3,9}*delx,graph_var{4,9}*delx];
maxv2=max(max(sizev2));
maxrange2=maxv2+(100-mod(maxv2,100));
plot(graph_var{1,9}*delx,'LineWidth',2,'Color',colorj(closestIndex1,:));
hold on
plot(graph_var{2,9}*delx,'LineWidth',2,'Color',colorj(closestIndex2,:));
hold on
plot(base_varRR{1,9}*delx,'LineWidth',2,'Color',colorj(closestIndex3,:));
hold on
plot(graph_var{3,9}*delx,'LineWidth',2,'Color',colorj(closestIndex4,:));
hold on
plot(graph_var{4,9}*delx,'LineWidth',2,'Color',colorj(closestIndex5,:));
ylabel({'Cell Area (μm^2)'},'FontSize',fa);
xlabel('Time (MCS)','FontSize',fa);
xticks([0 400 800 1200])
xticklabels({'0','400','800','1200'})
yticks([300 300+(maxrange1-300)/4 300+2*(maxrange1-300)/4  300+3*(maxrange1-300)/4 maxrange1])
yticklabels({num2str(300) num2str(300+(maxrange1-300)/4) num2str(300+2*(maxrange1-300)/4) num2str(300+3*(maxrange1-300)/4) num2str(maxrange1)})
ylim([300 maxrange1])

a = get(gca,'XTick');
set(gca,'XTick',a,'fontsize',fa-2)
b = get(gca,'YTick');
set(gca,'YTick',b,'fontsize',fa-2)
set(gca,'box','off') 
%title('TGF-β Treated','FontSize',fc,'FontWeight','bold');
leg=legend(num2str(graph_var{1,12}),num2str(graph_var{2,12}),baseval,num2str(graph_var{3,12}),num2str(graph_var{4,12}),'FontSize',fa-6, 'Interpreter', 'tex');
title(leg,legtitle);

% 6- 5 TGFB Treated Confluence
subplot(3,3,6)
plot(graph_var{1,10},'LineWidth',2,'Color',colorj(closestIndex1,:));
hold on
plot(graph_var{2,10},'LineWidth',2,'Color',colorj(closestIndex2,:));
hold on
plot(base_varRR{1,10},'LineWidth',2,'Color',colorj(closestIndex3,:));
hold on
plot(graph_var{3,10},'LineWidth',2,'Color',colorj(closestIndex4,:));
hold on
plot(graph_var{4,10},'LineWidth',2,'Color',colorj(closestIndex5,:));
xlabel('Time (MCS)','FontSize',fa);
xticks([0 400 800 1200])
xticklabels({'0','400','800','1200'})
yticks([0 0.25 0.5 0.75 1])
yticklabels({'0','0.25','0.5','0.75','1'})
ylabel({'Grid','Confluence'},'FontSize',fa);
a = get(gca,'XTick');
set(gca,'XTick',a,'fontsize',fa-2)
b = get(gca,'YTick');
set(gca,'YTick',b,'fontsize',fa-2)
set(gca,'box','off') 
%title('Confluence','FontSize',fb);

caxis([colormin colormax])
colormap(colorj)


%legend(graph_var{1,11},graph_var{2,11},basename,graph_var{3,11},graph_var{4,11});
%sgtitle('Control','FontSize',fc,'FontWeight','bold');

% 7-Steady State Pheno *****************************

phenoend=zeros(2,5);
title1=legtitle; %Specify Variable Title

xval=[graph_var{1,12},graph_var{2,12},1,graph_var{3,12},graph_var{4,12}];
phenoend(1,1)=graph_var{1,3}(end);
phenoend(1,2)=graph_var{2,3}(end);
phenoend(1,3)=base_varRR{1,3}(end);
phenoend(1,4)=graph_var{3,3}(end);
phenoend(1,5)=graph_var{4,3}(end);
phenoend(2,1)=graph_var{1,8}(end);
phenoend(2,2)=graph_var{2,8}(end);
phenoend(2,3)=base_varRR{1,8}(end);
phenoend(2,4)=graph_var{3,8}(end);
phenoend(2,5)=graph_var{4,8}(end);

subplot(3,3,7) % Specify plot number
plot(xval,phenoend(1,:),'-ob','MarkerFaceColor','b','LineWidth',2);
hold on
plot(xval,phenoend(2,:),'-or','MarkerFaceColor','r','LineWidth',2);
title('Steady State Phenotype','FontSize',fa); %'FontWeight','bold'
xlabel('Scaling Factor','FontSize',fa);
ylabel({'Cell','Phenotype'},'FontSize',fa);
yticks([0 0.25 0.5 0.75 1])
yticklabels({'0','0.25','0.5','0.75','1'});
xticks(xval)
xticklabels({num2str(xval(1)),num2str(xval(2)),num2str(xval(3)),num2str(xval(4)),num2str(xval(5))})
a = get(gca,'XTick');
set(gca,'XTick',a,'fontsize',fa-2)
b = get(gca,'YTick');
set(gca,'YTick',b,'fontsize',fa-2)
set(gca,'box','off') 
leg=legend('Control','TGF-β Treated','FontSize',fa-6, 'Interpreter', 'tex');


% 8-Steady State Size *****************************
title1=legtitle;% Specify variable title
sizeend=zeros(2,5);


sizeend(1,1)=graph_var{1,4}(end)*delx;
sizeend(1,2)=graph_var{2,4}(end)*delx;
sizeend(1,3)=base_varRR{1,4}(end)*delx;
sizeend(1,4)=graph_var{3,4}(end)*delx;
sizeend(1,5)=graph_var{4,4}(end)*delx;
sizeend(2,1)=graph_var{1,9}(end)*delx;
sizeend(2,2)=graph_var{2,9}(end)*delx;
sizeend(2,3)=base_varRR{1,9}(end)*delx;
sizeend(2,4)=graph_var{3,9}(end)*delx;
sizeend(2,5)=graph_var{4,9}(end)*delx;
maxv2=max(max(sizeend));
maxrange2=maxv2+(100-mod(maxv2,100));
minv2=min(min(sizeend));
minrange2=minv2+(100-mod(minv2,100));
rangediff=maxrange2-minrange2;

subplot(3,3,8) % Specify plot number
plot(xval,sizeend(1,:),'-ob','MarkerFaceColor','b','LineWidth',2);
hold on
plot(xval,sizeend(2,:),'-or','MarkerFaceColor','r','LineWidth',2);
title('Steady State Cell Size','FontSize',fa); %'FontWeight','bold'
xlabel('Scaling Factor','FontSize',fa);
ylabel({'Cell Area (μm^2)'},'FontSize',fa);
xticks(xval)
xticklabels({num2str(xval(1)),num2str(xval(2)),num2str(xval(3)),num2str(xval(4)),num2str(xval(5))})
yticks([300 300+(maxrange1-300)/4 300+2*(maxrange1-300)/4  300+3*(maxrange1-300)/4 maxrange1])
yticklabels({num2str(300) num2str(300+(maxrange1-300)/4) num2str(300+2*(maxrange1-300)/4) num2str(300+3*(maxrange1-300)/4) num2str(maxrange1)})
ylim([300 maxrange1])

a = get(gca,'XTick');
set(gca,'XTick',a,'fontsize',fa-2)
b = get(gca,'YTick');
set(gca,'YTick',b,'fontsize',fa-2)
leg=legend('Control','TGF-β Treated','FontSize',fa-6, 'Interpreter', 'tex');
set(gca,'box','off') 

% 9- SS confluence to 75%
subplot(3,3,9)
grid75=zeros(2,5);
gridindex=cell(2,5);

gridindex{1,1}=find(graph_var{1,5}>=0.75);
grid75(1,1)=gridindex{1,1}(1);

gridindex{1,2}=find(graph_var{2,5}>=0.75);
grid75(1,2)=gridindex{1,2}(1);

gridindex{1,3}=find(base_varRR{1,5}>=0.75);
grid75(1,3)=gridindex{1,3}(1);

gridindex{1,4}=find(graph_var{3,5}>=0.75);
grid75(1,4)=gridindex{1,4}(1);

gridindex{1,5}=find(graph_var{4,5}>=0.75);
grid75(1,5)=gridindex{1,5}(1);


gridindex{2,1}=find(graph_var{1,10}>=0.75);
grid75(2,1)=gridindex{2,1}(1);

gridindex{2,2}=find(graph_var{2,10}>=0.75);
grid75(2,2)=gridindex{2,2}(1);

gridindex{2,3}=find(base_varRR{1,10}>=0.75);
grid75(2,3)=gridindex{2,3}(1);

gridindex{2,4}=find(graph_var{3,10}>=0.75);
grid75(2,4)=gridindex{2,4}(1);

gridindex{2,5}=find(graph_var{4,10}>=0.75);
grid75(2,5)=gridindex{2,5}(1);

plot(xval,grid75(1,:),'-ob','MarkerFaceColor','b','LineWidth',2);
hold on
plot(xval,grid75(2,:),'-or','MarkerFaceColor','r','LineWidth',2);
title('Time to Reach 75% Confluence','FontSize',fa); %'FontWeight','bold'
xlabel('Scaling Factor','FontSize',fa);
ylabel({'Time (MCS)'},'FontSize',fa);
xticks(xval)
a = get(gca,'XTick');
set(gca,'XTick',a,'fontsize',fa-2)
b = get(gca,'YTick');
set(gca,'YTick',b,'fontsize',fa-2)
set(gca,'box','off') 


%% Visualizing 0 1 2.5 and 5 tgfb
fa=22;
fb=fa+2;
fc=fb+2;
colorj=jet;
var1=tgfb_varRR;
colormin=0;
colormax=5;
colorrange=colormin:(colormax-colormin)/255:colormax;
[minValue,closestIndex1] = min(abs(colorrange-0));
[minValue,closestIndex2] = min(abs(colorrange-1));
[minValue,closestIndex3] = min(abs(colorrange-2.5));
[minValue,closestIndex4] = min(abs(colorrange-5));
figure
subplot(1,3,1)

plot(var1{1,3},'LineWidth',2,'Color',colorj(closestIndex1,:));
hold on

plot(var1{2,3},'LineWidth',2,'Color',colorj(closestIndex2,:));
hold on

plot(var1{3,3},'LineWidth',2,'Color',colorj(closestIndex3,:));
hold on
plot(var1{4,3},'LineWidth',2,'Color',colorj(closestIndex4,:));
ylim([0 1])

a = get(gca,'XTick');
set(gca,'XTick',a,'fontsize',fa-2)
b = get(gca,'YTick');
set(gca,'YTick',b,'fontsize',fa-2)
xlabel('Time (MCS)','FontSize',fa);
xticks([0 400 800 1200])
xticklabels({'0','400','800','1200'})
yticks([0 0.25 0.5 0.75 1])
yticklabels({'0','0.25','0.5','0.75','1'})
ylabel('Cell Phenotype','FontSize',fa);
title('Phenotype','FontSize',fb);


subplot(1,3,2)
plot(var1{1,4}*delx,'LineWidth',2,'Color',colorj(closestIndex1,:));
hold on
plot(var1{2,4}*delx,'LineWidth',2,'Color',colorj(closestIndex2,:));
hold on
plot(var1{3,4}*delx,'LineWidth',2,'Color',colorj(closestIndex3,:));
hold on
plot(var1{4,4}*delx,'LineWidth',2,'Color',colorj(closestIndex4,:));
sizev1=[var1{1,4}*delx,var1{2,4}*delx,var1{3,4}*delx,var1{4,4}*delx];
maxv1=max(max([sizev1]));
maxrange1=maxv1+(100-mod(maxv1,100));
yticks([300 300+(maxrange1-300)/4 300+2*(maxrange1-300)/4  300+3*(maxrange1-300)/4 maxrange1])
yticklabels({num2str(300) num2str(300+(maxrange1-300)/4) num2str(300+2*(maxrange1-300)/4) num2str(300+3*(maxrange1-300)/4) num2str(maxrange1)})
ylim([300 maxrange1])
ylabel('Cell Area (μm^2)','FontSize',fa);
xlabel('Time (MCS)','FontSize',fa);
xticks([0 400 800 1200])
xticklabels({'0','400','800','1200'})
a = get(gca,'XTick');
set(gca,'XTick',a,'fontsize',fa-2)
b = get(gca,'YTick');
set(gca,'YTick',b,'fontsize',fa-2)
title('Cell Size','FontSize',fb);

subplot(1,3,3)

plot(var1{1,5},'LineWidth',2,'Color',colorj(closestIndex1,:));
hold on
plot(var1{2,5},'LineWidth',2,'Color',colorj(closestIndex2,:));
hold on
plot(var1{3,5},'LineWidth',2,'Color',colorj(closestIndex3,:));
hold on
plot(var1{4,5},'LineWidth',2,'Color',colorj(closestIndex4,:));
xlabel('Time (MCS)','FontSize',fa);
xticks([0 400 800 1200])
xticklabels({'0','400','800','1200'})
yticks([0 0.25 0.5 0.75 1])
yticklabels({'0','0.25','0.5','0.75','1'})
ylabel('Grid Confluence','FontSize',fa);
a = get(gca,'XTick');
set(gca,'XTick',a,'fontsize',fa-2)
b = get(gca,'YTick');
set(gca,'YTick',b,'fontsize',fa-2)
title('Confluence','FontSize',fb);

leg=legend(num2str(var1{1,1}),num2str(var1{2,1}),num2str(var1{3,1}),num2str(var1{4,1}),'FontSize',fa-2);
title(leg,'TGF-β');

caxis([colormin colormax])
colormap(jet)

%sgtitle('Comparing various doses of TGFB')

%% Visualize ECM Snapshots 7/15 *******************************************

load('test_100p_1500mcs_0tgfb_0.3tmax_0.06spread_div_0.0625_Jhalf9574.3821_extraA0.96_extraBT0.6_ZEBs0.8_R200s1.2_ZEB1_SNAIL1_kDAE1_PMN_scf0.5.mat');
%%
tic
fulltime=1500;
solTvar=zeros(100*100,fulltime);
solEvar=zeros(100*100,fulltime);
asmEvar=zeros(100*100,fulltime);
ecmTvar=zeros(100*100,fulltime);
for MCS=1:fulltime

maskG=ctag(:,MCS);
stateV=statevars{MCS};
NRC=size(csize{MCS},1);
k0T_val=zeros(1,NRC);
k0T_val(:)=0.06;
kT_val=zeros(1,NRC);
kT_val(:)=1.2;
nx=100*100;
MR200_val=zeros(sqrt(nx),sqrt(nx));
R2Val=stateV(30000+(5*NRC)+1: 30000+(6*NRC));
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

solEv=stateV(1:10000);
totalEv=stateV(10001:20000);
totalTv=stateV(20001:30000);
kbT=10.656*0.6;

Ncad_val=zeros(sqrt(nx),sqrt(nx));
NVal=stateV(30000+(7*NRC)+1: 30000+(8*NRC));
for i=1:NRC
    Ncad_temp=NVal(i);
    Ncad_val(find(maskG==i))=Ncad_temp;
end
Ncad_val=reshape(Ncad_val,nx,1);

param_kdae=1.225;
fdAE=5; % Scaling Factor
kdAE2=(720/500/25)*param_kdae;
kdAE1=kdAE2/fdAE; % from ECM degradation graph
kdAE=(kdAE2-kdAE1).*Ncad_val/3.1515 + kdAE1; % Uses Ncad so calculated here
kdAE=kdAE*1;

ecmTv=(kbT.*(totalTv+totalEv)+kdAE - sqrt( (kbT.*(totalTv+totalEv)+kdAE).^2 - 4.*kbT.^2.*totalTv.*totalEv))./(2*kbT);
solTv=totalTv-ecmTv;
asmEv=totalEv-ecmTv;

% storing variables
solTvar(:,MCS)=solTv;
solEvar(:,MCS)=solEv;
asmEvar(:,MCS)=asmEv;
ecmTvar(:,MCS)=ecmTv;

end
toc
%% Plotting at certain MCS for 4 concentrations
figure
colormap(jet);
for i=1:5
MCSval=510+(i-1)*50;
subplot(4,5,i)
imagesc(rot90(reshape(solTvar(:,MCSval),100,100)));
caxis([min(min(solTvar)), max(max(solTvar))]);
colorbar
title(['solT at ' num2str(MCSval)]);
subplot(4,5,i+5)
imagesc(rot90(reshape(solEvar(:,MCSval),100,100)));
caxis([min(min(solEvar)), max(max(solEvar))]);
colorbar
title(['solE at' num2str(MCSval)]);
subplot(4,5,i+10)
imagesc(rot90(reshape(asmEvar(:,MCSval),100,100)));
caxis([min(min(asmEvar)), max(max(asmEvar))]);
colorbar
title(['asmE at' num2str(MCSval)]);
subplot(4,5,i+15)
imagesc(rot90(reshape(ecmTvar(:,MCSval),100,100)));
caxis([min(min(ecmTvar)), max(max(ecmTvar))]);
colorbar
title(['ecmT at' num2str(MCSval)]);
end
sgtitle('5 TGFB');
%% Plotting ecmT

figure
colormap(jet);
totalEvar=asmEvar+ecmTvar;
MCSrange=[1000,1100,1200,1300]; % 1 figure with 600,700,800,900 and 1 with 1000,1100,1200,1300
for i=1:4
    
    MCSval=MCSrange(i);
    
    subplot(1,4,i)
    imagesc(rot90(reshape(ecmTvar(:,MCSval),100,100)));
    caxis([min(min(ecmTvar)), max(max(ecmTvar))]);
colorbar
title(['ecmT at MCS ' num2str(MCSval)]);
end
%% scratch
fa=24;
fb=fa+2;
fc=fb+2;
delx=2.5*2.5;
figure
subplot(1,3,1)
plot(wound_var{1,3},'LineWidth',2);
hold on
plot(wound_var{2,3},'LineWidth',2);
hold on
plot(wound_var{3,3},'LineWidth',3);
hold on
plot(wound_var{4,3},'LineWidth',2);
xlabel('Time (MCS)','FontSize',fa);
ylabel('Cell Phenotype','FontSize',fa);
title('Phenotype','FontSize',fb);
ylim([0 1]);
xticks([0 50 100 150])
xticklabels({'0','50','100','150'})
yticks([0 0.25 0.5 0.75 1])
yticklabels({'0','0.25','0.5','0.75','1'})
a = get(gca,'XTick');
set(gca,'XTick',a,'fontsize',fa-2)
b = get(gca,'YTick');
set(gca,'YTick',b,'fontsize',fa-2)

sizev=[wound_var{1,4}*delx; wound_var{2,4}*delx; wound_var{3,4}*delx; wound_var{4,4}*delx];
subplot(1,3,2)
plot(wound_var{1,4}*delx,'LineWidth',2);
hold on
plot(wound_var{2,4}*delx,'LineWidth',2);
hold on
plot(wound_var{3,4}*delx,'LineWidth',3);
hold on
plot(wound_var{4,4}*delx,'LineWidth',2);
ylabel('Cell Area (μm^2)','FontSize',fa);
xlabel('Time (MCS)','FontSize',fa);
xticks([0 50 100 150])
xticklabels({'0','50','100','150'})
title('Cell Size','FontSize',fb);
ylim([min(min(sizev))-10 max(max(sizev))+10]);
a = get(gca,'XTick');
set(gca,'XTick',a,'fontsize',fa-2)
b = get(gca,'YTick');
set(gca,'YTick',b,'fontsize',fa-2)

subplot(1,3,3)
plot(wound_var{1,5},'LineWidth',2);
hold on
plot(wound_var{2,5},'LineWidth',2);
hold on
plot(wound_var{3,5},'LineWidth',3);
hold on
plot(wound_var{4,5},'LineWidth',2);
ylim([0 1]);
xticks([0 50 100 150])
xticklabels({'0','50','100','150'})
ylim([0.5 1])
%yticks([0 0.25 0.5 0.75 1])
%yticklabels({'0','0.25','0.5','0.75','1'})
xlabel('Time (MCS)','FontSize',fa);
ylabel('Grid Confluence','FontSize',fa);
title('Confluence','FontSize',fb);
a = get(gca,'XTick');
set(gca,'XTick',a,'fontsize',fa-2)
b = get(gca,'YTick');
set(gca,'YTick',b,'fontsize',fa-2)

leg=legend(wound_var{:,6},'FontSize',fa-2);
title(leg,'Cell and Grid','FontSize',fb);

%% Binning- grid and over time
tic
load('test_100p_1201mcs_5tgfb_0.3tmax_0.06spread_Jhalf9574.3821_extraA0.96_extraBT0.6_ZEBscale0.8_R200scale1.2_ZEB1_SNAIL1_kDAE.scale1_TgfbRR.mat') %_inelasticity_1000
%%
endt=1200;
figure
dim=100; %dim x dim matrix
tp=endt; %timepoint
bin=5; % bin size
binsize=dim/bin;
cells=size(csize{1,tp},1);
loc1=reshape(ctag(:,tp),dim,dim);
loc=rot90(loc1); % Rotates matrix to match video orientation
size(loc);
finalsize=csize{1,tp};
finalp=phenotype{1,tp};
imagesc(loc)
colorbar
%colormap('jet');
title(['Cell #s, MCS ' num2str(tp)]);
b1=zeros(dim,dim);
p1=zeros(dim,dim);
%Centroid stuff
for i=1:cells
    [f1,f2]=find(loc == i);
    hold on
    xc=round(mean(f2),0);
    yc=round(mean(f1),0);
    text(xc,yc,[num2str(i)],'Color','w')
    b1(xc,yc)=finalsize(i);
    p1(xc,yc)=finalp(i);
end
%set(gca,'YDir','normal'); % start original at bottom right
b1=b1';%transposing makes it match video
p1=p1';
%b1 matrix with cell size at centroid
%
% figure
% imagesc(b1);
% colorbar
% figure
% imagesc(p1);
% colorbar
% colormap('jet');
% s([0 1]);
% % set(gca,'YDir','normal'); % start original at bottom right

% Size Binning
b2=zeros(dim,dim);
b1(b1 == 0) = NaN;
bvec=zeros(1,bin^2);
count=1;
for x=1:bin
    for y=1:bin
        %disp([num2str((20*(x-1)+1)),':', num2str((20*(x))),',',num2str((20*(y-1)+1)),':',num2str((20*(y)))]);
        val=mean(mean(b1((binsize*(y-1)+1):(binsize*(y)),(binsize*(x-1)+1):(binsize*(x))), 'omitnan'),'omitnan');
        %disp([num2str(y), ',' num2str(x), ': ', num2str(val)]);
        b2((binsize*(y-1)+1):(binsize*(y)),(binsize*(x-1)+1):(binsize*(x)))=val;
        bvec(count)=val;
        count=count+1;
    end
end
%b2 is binned cell size matrix
figure
imagesc(b2);
colorbar
colormap('jet');
caxis([0 200]);
title('Binned Cell Size','FontSize',18);
set(gca,'YDir','normal'); % start original at bottom right

% Phenotype Binning
p2=zeros(dim,dim);
p1(p1 == 0) = NaN;
binsize=dim/bin;
pvec=zeros(1,bin^2);
count=1;
for x=1:bin
    for y=1:bin
        %disp([num2str((20*(x-1)+1)),':', num2str((20*(x))),',',num2str((20*(y-1)+1)),':',num2str((20*(y)))]);
        val=mean(mean(p1((binsize*(y-1)+1):(binsize*(y)),(binsize*(x-1)+1):(binsize*(x))), 'omitnan'),'omitnan');
        %disp([num2str(y), ',' num2str(x), ': ', num2str(val)]);
        p2((binsize*(y-1)+1):(binsize*(y)),(binsize*(x-1)+1):(binsize*(x)))=val;
        pvec(count)=val;
        count=count+1;
    end
end
%p2 is binned cell phenotype matrix
% figure
% imagesc(p2);
% colorbar
% colormap('jet');
% caxis([0 1]);
% title('Binned Cell Phenotype','FontSize',18);
%set(gca,'YDir','normal'); % start original at bottom right

%centroid as average of x coordinates and y coordinates
%find centroid of each cell, replace centroid value with cell size.
%so you'll have eg. 100x100 matrix of 0s with 272 cells. do same with phenotype
% Junctional Forces
ctag_mat=reshape(ctag(:,tp),dim,dim);
[xjunc_cent, yjunc_cent, jx_cent, jy_cent, ...
    xjunc_vox, yjunc_vox, jx_vox_norm, jy_vox_norm, jx_vox, jy_vox,...
    jx_mat, jy_mat, neighbors_mat, Aconn] ...
    = calc_junc_from_ctag(ctag_mat, dim, dim);
% figure
% Q1=quiver(xjunc_cent,yjunc_cent,jx_cent,jy_cent);
% set(gca,'YDir','reverse'); % start original at bottom right
%
% binning Junctional forces
j1=zeros(dim,dim);
jvec=zeros(1,bin^2);
count=1;
for x=1:bin
    for y=1:bin
        %disp([num2str((20*(x-1)+1)),':', num2str((20*(x))),',',num2str((20*(y-1)+1)),':',num2str((20*(y)))]);
        xv1=jx_cent(xjunc_cent>=binsize*(x-1)+1 & xjunc_cent<binsize*(x) & yjunc_cent>=binsize*(y-1)+1 & yjunc_cent<binsize*(y));
        yv1=jy_cent(xjunc_cent>=binsize*(x-1)+1 & xjunc_cent<binsize*(x) & yjunc_cent>=binsize*(y-1)+1 & yjunc_cent<binsize*(y));
        val=mean(sqrt(xv1.^2+yv1.^2));
        %disp([num2str(y), ',' num2str(x), ': ', num2str(val)]);
        j1((binsize*(y-1)+1):(binsize*(y)),(binsize*(x-1)+1):(binsize*(x)))=val;
        jvec(count)=val;
        count=count+1;
    end
end
figure
imagesc(j1);
colorbar
colormap('jet');
% caxis([0 1]);
%title('Binned Junctional Forces','FontSize',18);
ax = gca;
ax.FontSize = 24;
% Mt2=jvec;


% Time Graphs with hard coded bin region finder
tic
phentime=zeros(1000,4);
sizetime=zeros(1000,4);
jforcetime=zeros(1000,4);
dim=100; %dim x dim matrix
finalt=endt; %timepoint
bin=5; % bin size
binsize=dim/bin;
%start forloop here
for tp=1:finalt
    disp(['MCS: ' num2str(tp)])
cells=size(csize{1,tp},1);
loc1=reshape(ctag(:,tp),dim,dim);
loc=rot90(loc1); % Rotates matrix to match video orientation
size(loc);
finalsize=csize{1,tp};
finalp=phenotype{1,tp};
b1=zeros(dim,dim);
p1=zeros(dim,dim);
%Centroid mapped over dim x dim matrix
for i=1:cells
    [f1,f2]=find(loc == i);
    xc=round(mean(f2),0);
    yc=round(mean(f1),0);
    b1(xc,yc)=finalsize(i);
    p1(xc,yc)=finalp(i);
end
b1=b1';%transposing makes it match video - b1 matrix with cell size at centroid
p1=p1'; %p1 matrix for phenotype

ctag_mat=reshape(ctag(:,tp),dim,dim);
[xjunc_cent, yjunc_cent, jx_cent, jy_cent, ...
    xjunc_vox, yjunc_vox, jx_vox_norm, jy_vox_norm, jx_vox, jy_vox,...
    jx_mat, jy_mat, neighbors_mat, Aconn] ...
    = calc_junc_from_ctag(ctag_mat, dim, dim);
%xcent,ycent= locations, jxcent,jycent=forces
%Binning
b2=zeros(dim,dim);
b1(b1 == 0) = NaN;
p2=zeros(dim,dim);
p1(p1 == 0) = NaN;
j1=zeros(dim,dim);
binsize=dim/bin;
pvec=zeros(1,bin^2);
bvec=zeros(1,bin^2);
jvec=zeros(1,bin^2);
count=1;
for x=1:bin
    for y=1:bin
        %size
        val1=mean(mean(b1((binsize*(y-1)+1):(binsize*(y)),(binsize*(x-1)+1):(binsize*(x))), 'omitnan'),'omitnan');
        b2((binsize*(y-1)+1):(binsize*(y)),(binsize*(x-1)+1):(binsize*(x)))=val1;
        bvec(count)=val1;        
        %phenotype
        val2=mean(mean(p1((binsize*(y-1)+1):(binsize*(y)),(binsize*(x-1)+1):(binsize*(x))), 'omitnan'),'omitnan');
        p2((binsize*(y-1)+1):(binsize*(y)),(binsize*(x-1)+1):(binsize*(x)))=val2;
        pvec(count)=val2;
        %junctional forces  
        xv1=jx_cent(xjunc_cent>=binsize*(x-1)+1 & xjunc_cent<binsize*(x) & yjunc_cent>=binsize*(y-1)+1 & yjunc_cent<binsize*(y));
        yv1=jy_cent(xjunc_cent>=binsize*(x-1)+1 & xjunc_cent<binsize*(x) & yjunc_cent>=binsize*(y-1)+1 & yjunc_cent<binsize*(y));
        val3=mean(sqrt(xv1.^2+yv1.^2),'omitnan');
        j1((binsize*(y-1)+1):(binsize*(y)),(binsize*(x-1)+1):(binsize*(x)))=val3;
        jvec(count)=val3;
        count=count+1;
    end
end
%reshape to 5x5
pvec=reshape(pvec,5,5)';
bvec=reshape(bvec,5,5)';
jvec=reshape(jvec,5,5)';

pcenter=pvec(13);
bcenter=bvec(13);
jcenter=jvec(13);

pedge=mean([pvec(1,1),pvec(1,5),pvec(5,1),pvec(5,5)],'omitnan');
bedge=mean([bvec(1,1),bvec(1,5),bvec(5,1),bvec(5,5)],'omitnan');
jedge=mean([jvec(1,1),jvec(1,5),jvec(5,1),jvec(5,5)],'omitnan');

pmid=mean([pvec(2,2),pvec(2,3),pvec(2,4),pvec(3,2),pvec(3,4),pvec(4,2),pvec(4,3),pvec(4,4)],'omitnan');
bmid=mean([bvec(2,2),bvec(2,3),bvec(2,4),bvec(3,2),bvec(3,4),bvec(4,2),bvec(4,3),bvec(4,4)],'omitnan');
jmid=mean([jvec(2,2),jvec(2,3),jvec(2,4),jvec(3,2),jvec(3,4),jvec(4,2),jvec(4,3),jvec(4,4)],'omitnan');

pbor=mean([pvec(1,2),pvec(1,3),pvec(1,4),pvec(2,1),pvec(2,5),pvec(3,1),pvec(3,5),pvec(4,1),pvec(4,5),pvec(5,2),pvec(5,3),pvec(5,4)],'omitnan');
bbor=mean([bvec(1,2),bvec(1,3),bvec(1,4),bvec(2,1),bvec(2,5),bvec(3,1),bvec(3,5),bvec(4,1),bvec(4,5),bvec(5,2),bvec(5,3),bvec(5,4)],'omitnan');
jbor=mean([jvec(1,2),jvec(1,3),jvec(1,4),jvec(2,1),jvec(2,5),jvec(3,1),jvec(3,5),jvec(4,1),jvec(4,5),jvec(5,2),jvec(5,3),jvec(5,4)],'omitnan');

phentime(tp,1)=pcenter;
sizetime(tp,1)=bcenter;
jforcetime(tp,1)=jcenter;

phentime(tp,2)=pmid;
sizetime(tp,2)=bmid;
jforcetime(tp,2)=jmid;

phentime(tp,3)=pbor;
sizetime(tp,3)=bbor;
jforcetime(tp,3)=jbor;

phentime(tp,4)=pedge;
sizetime(tp,4)=bedge;
jforcetime(tp,4)=jedge;

end
%%
fa=23;
figure %Center
 subplot(1,3,1)
plot(phentime(:,1),'LineWidth',3,'Color','r');
hold on
plot(phentime(:,2),'LineWidth',3,'Color','#EDB120');
hold on
plot(phentime(:,3),'LineWidth',3,'Color','#4DBEEE');
hold on
plot(phentime(:,4),'LineWidth',3,'Color','#77AC30');
title('Cell Phenotype','FontSize',24);
xlim([0,endt])
xticks([0 400 800 1200])
xticklabels({'0','400','800','1200'})
ylim([0,1])
%legend('Center','Middle','Border','Corner (Edges)','FontSize',30);
xlabel('Time (MCS)','FontSize',fa);
ylabel('Phenotype','FontSize',fa);
a = get(gca,'XTick');
set(gca,'XTick',a,'fontsize',fa-2)
b = get(gca,'YTick');
set(gca,'YTick',b,'fontsize',fa-2)
set(gca,'box','off') 
%set(gca,'box','off','tickdir','out','LineWidth',2,'ycolor','w')
subplot(1,3,2)
plot(sizetime(:,1)*delx,'LineWidth',3,'Color','r');
hold on
plot(sizetime(:,2)*delx,'LineWidth',3,'Color','#EDB120');
hold on
plot(sizetime(:,3)*delx,'LineWidth',3,'Color','#4DBEEE');
hold on
plot(sizetime(:,4)*delx,'LineWidth',3,'Color','#77AC30');
maxv1=max(max([sizetime*delx]));
maxrange1=maxv1+(100-mod(maxv1,100));
yticks([300 300+(maxrange1-300)/4 300+2*(maxrange1-300)/4  300+3*(maxrange1-300)/4 maxrange1])
yticklabels({num2str(300) num2str(300+(maxrange1-300)/4) num2str(300+2*(maxrange1-300)/4) num2str(300+3*(maxrange1-300)/4) num2str(maxrange1)})
ylim([300 maxrange1])
title('Cell Size','FontSize',24)
ylabel('Cell Area (μm^2)','FontSize',fa);
xlim([0,endt])
%ylim([0,200])
xlabel('Time (MCS)','FontSize',fa);
xticks([0 400 800 1200])
xticklabels({'0','400','800','1200'})
a = get(gca,'XTick');
set(gca,'XTick',a,'fontsize',fa-2)
b = get(gca,'YTick');
set(gca,'YTick',b,'fontsize',fa-2)
set(gca,'box','off') 

subplot(1,3,3)
plot(jforcetime(:,1),'LineWidth',3,'Color','r');
hold on
plot(jforcetime(:,2),'LineWidth',3,'Color','#EDB120');
hold on
plot(jforcetime(:,3),'LineWidth',3,'Color','#4DBEEE');
hold on
plot(jforcetime(:,4),'LineWidth',3,'Color','#77AC30');
title('Junctional Forces','FontSize',24)
xlim([0,endt])
%ylim([0,50000])
legend('Center','Middle','Border','Corner (Edges)');
%legend('Center','Middle','Border','Corner (Edges)','FontSize',30);
xlabel('Time (MCS)','FontSize',fa);
xticks([0 400 800 1200])
xticklabels({'0','400','800','1200'})
a = get(gca,'XTick');
set(gca,'XTick',a,'fontsize',fa-2)
b = get(gca,'YTick');
set(gca,'YTick',b,'fontsize',fa-2)
set(gca,'box','off') 
max(max(jforcetime))

toc


%% visualizing grid over time

var1=sizeAvg0;
maxv=max(max(var1));
%maxv=1;
for i=1:12
    figure
    
    imagesc(reshape(var1(:,i),dim,dim))
    title(['Average Size for param: ' parameter ' at MCS: ' num2str(i*dim)]);
    colorbar
    colormap('jet');
    caxis([0 maxv]);
    pause(0.5)
end
%%
function files =GetNames(filetype)
files2 = dir(fullfile(['**/*'], filetype));
vec=[];
for i=1:size(files2,1)
    vec=[vec ',' [files2(i).name]];
end
%vec2=convertCharsToStrings(vec);
files=strsplit(vec,',');
files=files(2:end);
end