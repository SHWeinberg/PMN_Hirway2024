load('test_100p_3001mcs_0tgfb_0.06spread0div_ECMscale1__21EMTresPer_METproneCent_BaseGrid.mat');
%%
tic
endt=size(csize,2);
binsizemat=cell(endt,1);
degree_mat=cell(endt,1);
bin_mat=cell(endt,1);
cells=size(csize{1},1);

for i=1:endt
    % i=650;
    G=graph(conn_mat{i});
    % plot(G)
    % x0=10;
    % y0=150;
    % width=600;
    % height=500;
    % set(gcf,'position',[x0,y0,width,height])
    % title(['MCS: ' num2str(i)]);
    D = degree(G);
    degree_mat{i}=D;
    [bins,binsizes] = conncomp(G);
    bin_mat{i}=bins;
    binsizemat{i}=binsizes;
    %pause(0.1)
end
% %%
% figure
% i=465;
% G=graph(conn_mat{i});
% plot(G)
% title(['MCS: ' num2str(i) ', bins: ' num2str(binsizemat{i})]);

%
cellbin=zeros(endt,cells);
cellnum=1;
for j=1:cells
    for i=1:endt
        cellbin(i,j)=bin_mat{i}(j);
    end
end
%
% figure
% i=95;
% plot(cellbin(:,i));
% hold on
% plot(cellbin(:,i+1));
% hold on
% plot(cellbin(:,i+2));
% hold on
% plot(cellbin(:,i+3));
% hold on
% plot(cellbin(:,i+4));
% title(['Cells ' num2str(i) ' to ' num2str(i+5)]);
% legend(num2str(i),num2str(i+1),num2str(i+2),num2str(i+3),num2str(i+4));

%
phenotypemat=zeros(endt,cells);
for j=1:cells
    for i=1:endt
        phenotypemat(i,j)=phenotype{i}(j);
    end
end
% Average periphery phenotype over time
peripheryP=phenotypemat(:,1:92);
peripheryP=mean(peripheryP');
centerP=phenotypemat(:,93:100);
centerP=mean(centerP');
figure
plot(centerP);
hold on
plot(peripheryP);
title('Grid Phenotype');
xlabel('MCS');
ylabel('Phenotype');
ylim([-0.1 1.1]);
xlim([0 endt]);
legend('Center','Periphery');
%
peribin=zeros(endt,92);
for i=1:endt
    peribin(i,:)=bin_mat{i}(1:92);
end
centerbin=zeros(endt,8);
for i=1:endt
    centerbin(i,:)=bin_mat{i}(93:100);
end
N2=zeros(endt,max(max(peribin)));

% Track invasion

% for i=1:endt
%     for n=1:max(max(peribin))
%         N2(i,n)=sum(peribin(i,:)==n);
%     end
% end
periphery=zeros(endt,1);
clustr=zeros(endt,1);
for i=1:endt
    clustr(i)=size(centerbin,2);
    periphery(i)=size(peribin,2);
    for c=1:size(centerbin,2)
        
    if(sum((peribin(i,:)==centerbin(i,c)))>0)
        clustr(i)=clustr(i)-1;
        periphery(i)=periphery(i)+1;
    end
    end
end

figure 
yyaxis left
plot(clustr)
ylim([0,10]);
ylabel('Number of Cells in Center');
yyaxis right
plot(periphery)
ylim([92,102]);
xlim([0 endt]);
ylabel('Number of Cells in Periphery');
title('Center and Periphery Size');
legend('Center','Periphery');
xlabel('MCS');
toc

%%
bin1=zeros(endt,1);
bin2=zeros(endt,1);
for i=1:endt
    bin1(i)=binsizemat{i}(1);
    if size(binsizemat{i},2)>1
        bin2(i)=binsizemat{i}(2);
    else
        bin2(i)=0;
    end
end
figure
plot(bin1);
hold on
plot(bin2);
legend('bin1','bin2');


%% Track mesenchymal center
centerbox=zeros(1681,endt);
uniq1=cell(endt,1);
uniq2=zeros(endt,1);
for i=1:endt
    reshaped1=reshape(ctag(:,i),100,100);
    reshaped2=reshape(reshaped1(30:70,30:70),1681,1);
centerbox(:,i)=reshaped2;
uniq1{i}=unique(centerbox(:,i));
uniq2(i)=size(find(uniq1{i}>92),1);
end

figure
plot(uniq2);
xlabel('Time (MCS)');
ylabel('Center Size');

% N=zeros(endt,max(max(centerbin)));
% 
% for i=1:endt
%     for n=1:max(max(centerbin))
%         N(i,n)=sum(centerbin(i,:)==n);
%     end
% end
% 
% figure
% for i=1:max(max(centerbin))
% plot(N(:,i))
% hold on
% end
% ylim([0,10]);
%legend('Bin 1', 'Bin 2','Bin 3','Bin 4','Bin 5','Bin 6','Bin 7','Bin 8');
