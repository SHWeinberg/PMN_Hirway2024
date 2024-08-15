%1/13- SUH
tic
load('test_100p_500mcs_0tgfb_0.3tmax_0.06spread_targetscale2_inelasticity_1000_Jhalf11489.2586_dividing1.25_extra_A_BT0.93_ZEBscale0.8_R200scale1.2_.mat') %_inelasticity_1000
endt=3000;
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
%
figure %Center
 %subplot(1,3,1)
plot(phentime(:,1),'LineWidth',3,'Color','r');
hold on
plot(phentime(:,2),'LineWidth',3,'Color','#EDB120');
hold on
plot(phentime(:,3),'LineWidth',3,'Color','#4DBEEE');
hold on
plot(phentime(:,4),'LineWidth',3,'Color','#77AC30');
title('Cell Phenotype over time','FontSize',18);
xlim([0,endt])
ylim([0,1])
%legend('Center','Middle','Border','Corner (Edges)','FontSize',30);
xlabel('Time (MCS)','FontSize',18);
ax = gca;
ax.FontSize = 18;
%set(gca,'box','off','tickdir','out','LineWidth',2,'ycolor','w')
figure %subplot(1,3,2)
plot(sizetime(:,1),'LineWidth',3,'Color','r');
hold on
plot(sizetime(:,2),'LineWidth',3,'Color','#EDB120');
hold on
plot(sizetime(:,3),'LineWidth',3,'Color','#4DBEEE');
hold on
plot(sizetime(:,4),'LineWidth',3,'Color','#77AC30');
title('Cell Size over time','FontSize',18)
xlim([0,endt])
%ylim([0,200])
legend('Center','Middle','Border','Corner (Edges)');
xlabel('Time (MCS)');

figure %subplot(1,3,3)
plot(jforcetime(:,1),'LineWidth',3,'Color','r');
hold on
plot(jforcetime(:,2),'LineWidth',3,'Color','#EDB120');
hold on
plot(jforcetime(:,3),'LineWidth',3,'Color','#4DBEEE');
hold on
plot(jforcetime(:,4),'LineWidth',3,'Color','#77AC30');
title('Junctional Forces over time','FontSize',18)
xlim([0,endt])
%ylim([0,50000])
%legend('Center','Middle','Border','Corner (Edges)','FontSize',30);
xlabel('Time (MCS)','FontSize',18);
%set(gca,'box','off','tickdir','out','LineWidth',2,'ycolor','w')
ax = gca;
ax.FontSize = 18;
max(max(jforcetime))

toc


%% Gettting the U shape curve
%%
close all
al=0.05; % .05
n=2; %2
Jmin=4819; %4819
Jhalf0=Jmin/((al^-1-1)^(1/n));
Tmax=10;
tmax0=1;
Jmin= Jhalf0*(al^-1 -1) ^(1/n) / (2*Tmax/tmax0 - 1)^(1/n)
Jhalf=Jmin/((al^-1-1)^(1/n));
% Junctional Force analysis
Ejvec=reshape(Et2,5,5);
Pjvec=reshape(Pt2,5,5);
Mjvec=reshape(Mt2,5,5);
%Tmax=2; %1
Js=1.6;
Fmax=Tmax/Js; %1.6
al=0.05; % .05
n=2; %2
%Jmin= 2.7823e+03; %4819
Jhalf=Jmin/((al^-1-1)^(1/n));
EFmat=zeros(5,5);
PFmat=zeros(5,5);
MFmat=zeros(5,5);
for i=1:25
EFmat(i)= Fmax/(1+(Ejvec(i)/Jhalf)^2);
PFmat(i)= Fmax/(1+(Pjvec(i)/Jhalf)^2);
MFmat(i)= Fmax/(1+(Mjvec(i)/Jhalf)^2);
end
%
figure
plot(linspace(1,5,5),EFmat(3,:)*Js,'LineWidth',2)
hold on
plot(linspace(1,5,5),PFmat(3,:)*Js,'LineWidth',2)
hold on
plot(linspace(1,5,5),MFmat(3,:)*Js,'LineWidth',2)
legend('E','P','M');
title(['Cross-Section, Jmin:' num2str(Jmin) ', Tmax:' num2str(Tmax), ', targetscalex2, inelasticity1000'])
%ylim([0,0.2]);

%% *******************************************************************
% 11/22/2020
% Quantifying Extracellular EMT CPM
filename='test_100p_304mcs_5tgfb_0.25tmax_0.06spread_targetscale2_inelasticity_1000_Jhalf10212.6743_dividing1.25_extra_A_BT0.963_ZEBscale0.8_R200scale1.2_.mat';
load(filename);
%% 
tic
extrafilename=['extracellular_' filename];
endv=size(statevars,2);
myVideo = VideoWriter(extrafilename); %open video file
myVideo.FrameRate = 15;  %can adjust this, 5 - 10 works well for me
open(myVideo)
figure('Renderer', 'painters', 'Position', [5 5 900 700])
solTv=[];
solEv=[];
asmEv=[];
ecmTv=[];
Ncad1=[];
maxsolT=[];
maxsolE=[];
maxasmE=[];
maxecmT=[];
for time=1:endv
    a=statevars{1,time};
    maxsolT=[maxsolT; a(00001:10000)];
    maxsolE=[maxsolE; a(10001:20000)];
    maxasmE=[maxasmE; a(20001:30000)];
    maxecmT=[maxecmT; a(30001:40000)];
end

for i=1:endv
    subplot(2,2,1)
    a=statevars{1,i};
    b=a(00001:10000);
    c=reshape(b,100,100);
    imagesc(c);
    title(['solT at MCS:' num2str(i)]);
    colorbar
    caxis([0 max(maxsolT)]);
    
    subplot(2,2,2)
    a=statevars{1,i};
    b=a(10001:20000);
    c=reshape(b,100,100);
    imagesc(c);
    title('solE');
    colorbar
    caxis([0 max(maxsolE)]);
    
    subplot(2,2,3)
    a=statevars{1,i};
    b=a(20001:30000);
    c=reshape(b,100,100);
    imagesc(c);
    title('asmE');
    colorbar
    caxis([0 max(maxasmE)]);
    
    subplot(2,2,4)
    a=statevars{1,i};
    b=a(30001:40000);
    c=reshape(b,100,100);
    imagesc(c);
    title('ecmT');
    colorbar
    caxis([0 max(maxecmT)]);
    
    
    pause(0.001);
    colormap(jet);
    a=statevars{1,i};
    solTv=[solTv,a(345)];
    solEv=[solEv,a(10345)];
    asmEv=[asmEv,a(20345)];
    ecmTv=[ecmTv,a(30345)];
    Ncad1=[Ncad1,a(end)];
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
end
close(myVideo);
toc
%%
t=[1:304];
figure
subplot(3,2,1)
plot(t,solTv);
title('solT');

subplot(3,2,2)
plot(t,solEv);
title('solE');

subplot(3,2,3)
plot(t,asmEv);
title('asmE');

subplot(3,2,4)
plot(t,ecmTv);
title('ecmT');

subplot(3,2,5)
plot(t,Ncad1);
title('Ncad');