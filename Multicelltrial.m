%% Initial given parameters of cell in E state
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
init.exo=0; % Exogenous TGFB in the system
init.Tbas=0.06;
init.mir2=0.0002;

init.mask=1;
init.time=[0:1:24*20];% hours

%%
n =5; %  dimension  E cells
b=1;  % outer most boundary of no cells
m=1; % middle layer of M cells
dim=n+2*b+2*m; % Dimension of full matrix

tic
disp(['Running Cell Simulation'])
f1=Multicell(init,n,b,m);
toc



    
   
%% Heat Map of TGFB
f2=f1(:,1:dim^2);
tspan=init.time;
xv=linspace(1,dim,dim);
yv=linspace(1,dim,dim);
maxv=max(max(f2));
filename=['Diffusion TGFB (' num2str(init.exo) ' exo).gif'];
figure(1)
for i=1:size(tspan,2)
v=reshape(f2(i,:),dim,dim);
colormap('jet')
imagesc(xv,yv,v)
title(['TGFB Diffusion at time: ' num2str(i-1) ' hours'])
colorbar
caxis([0 maxv])
drawnow
frame=getframe(1);
im=frame2im(frame);
[imind,cm]=rgb2ind(im,256);
if i==1
    imwrite(imind,cm,filename,'gif','Loopcount',inf);
else
    imwrite(imind,cm,filename,'gif','WriteMode','append');
end

pause(0.01) %in seconds 

end

%% Time Plots of TGFB and of Ncad
tspan=init.time;
f2=f1(:,1:dim^2);
x=tspan;
m1=(dim*(b+m))+b+m;
c1=round(dim^2/2);
yc=f2(:,c1); % Center Cell
ym=f2(:,m1); %Mesenchymal
yb=f2(:,1); %Boundary side
figure
plot(x,yc,'g');
hold on
plot(x,ym,'r');
hold on 
plot(x,yb,'b');
title('TGFB Concentration over Time');
legend('Cell Center','MesenchymalCorner','Boundary');
xlabel('Hours')
ylabel('uM')

% Time Plots of N cad
tspan=init.time;
f2=f1(:,end-dim^2+1:end);
x=tspan;
m1=(dim*(b+m))+b+m;
c1=round(dim^2/2);
yc=f2(:,c1); % Center Cell
ym=f2(:,m1); %Mesenchymal
yb=f2(:,1); %Boundary side
figure
plot(x,yc,'g');
hold on
plot(x,ym,'r');
hold on 
plot(x,yb,'b');
title('N cad Concentration over Time');
legend('Cell Center','Mesenchymal Corner','Boundary');
xlabel('Hours')
ylabel('uM')
%% Heat Map of N cad
f2=f1(:,end-dim^2+1:end);
tspan=init.time;
xv=linspace(1,dim,dim);
yv=linspace(1,dim,dim);
maxv=max(max(f2));
filename=['Diffusion Ncad (' num2str(init.exo) ' exo).gif'];
figure(1)
for i=1:size(tspan,2)
v=reshape(f2(i,:),dim,dim);
colormap('jet')
imagesc(xv,yv,v)
title(['Ncad at Time: ' num2str(i-1) ' hours'])
colorbar
caxis([0 maxv])
drawnow
frame=getframe(1);
im=frame2im(frame);
[imind,cm]=rgb2ind(im,256);
if i==1
    imwrite(imind,cm,filename,'gif','Loopcount',inf);
else
    imwrite(imind,cm,filename,'gif','WriteMode','append');
end

pause(0.001) %in seconds 

end
%% Testing Case 1 Large epithelial (5E 1M)- No inhibitor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

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
init.exo=1; % Exogenous TGFB in the system
init.Tbas=0.06;
init.mir2=0.0002;

init.mask=1;
init.time=[0:1:24*20];% hours


%
n =5; %  dimension  E cells
b=1;  % outer most boundary of no cells
m=1; % middle layer of M cells
dim=n+2*b+2*m; % Dimension of full matrix

tic
disp(['Running Cell Simulation'])
f1=Multicell(init,n,b,m);
toc

% Time Plots of TGFB and of Ncad
tspan=init.time;
f2=f1(:,1:dim^2);
x=tspan;
m1=(dim*(b+m))+b+m;
c1=round(dim^2/2);
yc=f2(:,c1); % Center Cell
ym=f2(:,m1); %Mesenchymal
yb=f2(:,1); %Boundary side
figure
plot(x,yc,'g');
hold on
plot(x,ym,'r');
hold on 
plot(x,yb,'b');
title('TGFB Concentration over Time');
legend('Cell Center','MesenchymalCorner','Boundary');
xlabel('Hours')
ylabel('uM')
movegui('west');
% Time Plots of N cad
tspan=init.time;
f2=f1(:,end-dim^2+1:end);
x=tspan;
m1=(dim*(b+m))+b+m;
c1=round(dim^2/2);
yc=f2(:,c1); % Center Cell
ym=f2(:,m1); %Mesenchymal
yb=f2(:,1); %Boundary side
figure
plot(x,yc,'g');
hold on
plot(x,ym,'r');
hold on 
plot(x,yb,'b');
title('N cad Concentration over Time');
legend('Cell Center','Mesenchymal Corner','Boundary');
xlabel('Hours')
ylabel('uM')
movegui('east');

%% Testing Case 2 Large Mesenchymal (1E 5M)- No inhibitor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%

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
init.exo=1; % Exogenous TGFB in the system
init.Tbas=0.06;
init.mir2=0.0002; %0.0008460221

init.mask=1;
init.time=[0:1:24*20];% hours

%
n =1; %  dimension  E cells
b=1;  % outer most boundary of no cells
m=5; % middle layer of M cells
dim=n+2*b+2*m; % Dimension of full matrix

tic
disp(['Running Cell Simulation'])
f1=Multicell(init,n,b,m);
toc

% Time Plots of TGFB and of Ncad
tspan=init.time;
f2=f1(:,1:dim^2);
x=tspan;
m1=(dim*(b+m))+b+m;
c1=round(dim^2/2);
yc=f2(:,c1); % Center Cell
ym=f2(:,m1); %Mesenchymal
yb=f2(:,1); %Boundary side
figure
plot(x,yc,'g');
hold on
plot(x,ym,'r');
hold on 
plot(x,yb,'b');
title('TGFB Concentration over Time');
legend('Cell Center','MesenchymalCorner','Boundary');
xlabel('Hours')
ylabel('uM')
movegui('west');
% Time Plots of N cad
tspan=init.time;
f2=f1(:,end-dim^2+1:end);
x=tspan;
m1=(dim*(b+m))+b+m;
c1=round(dim^2/2);
yc=f2(:,c1); % Center Cell
ym=f2(:,m1); %Mesenchymal
yb=f2(:,1); %Boundary side
figure
plot(x,yc,'g');
hold on
plot(x,ym,'r');
hold on 
plot(x,yb,'b');
title('N cad Concentration over Time');
legend('Cell Center','Mesenchymal Corner','Boundary');
xlabel('Hours')
ylabel('uM')
movegui('east');


%% Testing Case 3 Large epithelial (5E 1M)- with inhibitor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

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
init.exo=1; % Exogenous TGFB in the system
init.Tbas=0.06;
init.mir2=0.0008460221;

init.mask=1;
init.time=[0:1:24*20];% hours


%
n =5; %  dimension  E cells
b=1;  % outer most boundary of no cells
m=1; % middle layer of M cells
dim=n+2*b+2*m; % Dimension of full matrix

tic
disp(['Running Cell Simulation'])
f1=Multicell(init,n,b,m);
toc

% Time Plots of TGFB and of Ncad
tspan=init.time;
f2=f1(:,1:dim^2);
x=tspan;
m1=(dim*(b+m))+b+m;
c1=round(dim^2/2);
yc=f2(:,c1); % Center Cell
ym=f2(:,m1); %Mesenchymal
yb=f2(:,1); %Boundary side
figure
plot(x,yc,'g');
hold on
plot(x,ym,'r');
hold on 
plot(x,yb,'b');
title('TGFB Concentration over Time');
legend('Cell Center','MesenchymalCorner','Boundary');
xlabel('Hours')
ylabel('uM')
movegui('west');
% Time Plots of N cad
tspan=init.time;
f2=f1(:,end-dim^2+1:end);
x=tspan;
m1=(dim*(b+m))+b+m;
c1=round(dim^2/2);
yc=f2(:,c1); % Center Cell
ym=f2(:,m1); %Mesenchymal
yb=f2(:,1); %Boundary side
figure
plot(x,yc,'g');
hold on
plot(x,ym,'r');
hold on 
plot(x,yb,'b');
title('N cad Concentration over Time');
legend('Cell Center','Mesenchymal Corner','Boundary');
xlabel('Hours')
ylabel('uM')
movegui('east');

%% Testing Case 4 Large Mesenchymal (1E 5M)- with inhibitor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%

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
init.exo=1; % Exogenous TGFB in the system
init.Tbas=0.06;
init.mir2=0.0008460221;

init.mask=1;
init.time=[0:1:24*20];% hours

%
n =1; %  dimension  E cells
b=1;  % outer most boundary of no cells
m=5; % middle layer of M cells
dim=n+2*b+2*m; % Dimension of full matrix

tic
disp(['Running Cell Simulation'])
f1=Multicell(init,n,b,m);
toc

% Time Plots of TGFB and of Ncad
tspan=init.time;
f2=f1(:,1:dim^2);
x=tspan;
m1=(dim*(b+m))+b+m;
c1=round(dim^2/2);
yc=f2(:,c1); % Center Cell
ym=f2(:,m1); %Mesenchymal
yb=f2(:,1); %Boundary side
figure
plot(x,yc,'g');
hold on
plot(x,ym,'r');
hold on 
plot(x,yb,'b');
title('TGFB Concentration over Time');
legend('Cell Center','MesenchymalCorner','Boundary');
xlabel('Hours')
ylabel('uM')
movegui('west');
% Time Plots of N cad
tspan=init.time;
f2=f1(:,end-dim^2+1:end);
x=tspan;
m1=(dim*(b+m))+b+m;
c1=round(dim^2/2);
yc=f2(:,c1); % Center Cell
ym=f2(:,m1); %Mesenchymal
yb=f2(:,1); %Boundary side
figure
plot(x,yc,'g');
hold on
plot(x,ym,'r');
hold on 
plot(x,yb,'b');
title('N cad Concentration over Time');
legend('Cell Center','Mesenchymal Corner','Boundary');
xlabel('Hours')
ylabel('uM')
movegui('east');




%% NCad vs TGFB Case 1
% Testing Case 1 Large epithelial (5E 1M)- No inhibitor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

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

init.mask=1;
init.time=[0:1:24*20];% hours

%
n =5; %  dimension  E cells
b=1;  % outer most boundary of no cells
m=1; % middle layer of M cells
dim=n+2*b+2*m; % Dimension of full matrix

exoc=[];
stored=[];
for i=1:81
init.exo=((i-1)/32); %Exogenous concentration of TGF-? 0?5 unit ******Can be changed to 3.0
exoc=[exoc; init.exo];
tic
disp(['Running Cell Simulation at iter: ' num2str(i)])
f1=Multicell(init,n,b,m);
toc
s1=f1(end,end-dim^2+1:end);
 stored=[stored;s1];
end
m1=(dim*(b+m))+b+m;
c1=round(dim^2/2);
yc1=stored(:,c1); % Center Cell
ym1=stored(:,m1); %Mesenchymal
yb1=stored(:,1); %Boundary side

figure
plot(exoc,yc);
hold on
plot(exoc,ym);
hold on
plot(exoc,yb);
ylabel('N-cadherin')
xlabel('TGFB')
legend('Cell Center','Mesenchymal Corner','Boundary','Location','southeast')

%% NCad vs TGFB Case 2
% Testing Case 2 Large epithelial (1E 5M)- No inhibitor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

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

init.mask=1;
init.time=[0:1:24*20];% hours


%
n =1; %  dimension  E cells
b=1;  % outer most boundary of no cells
m=5; % middle layer of M cells
dim=n+2*b+2*m; % Dimension of full matrix

exoc=[];
stored=[];
for i=1:81
init.exo=((i-1)/20); %Exogenous concentration of TGF-? 0?5 unit ******Can be changed to 3.0
exoc=[exoc; init.exo];
tic
disp(['Running Cell Simulation with exo: ' num2str(init.exo)])
f1=Multicell(init,n,b,m);
toc
s1=f1(end,end-dim^2+1:end);
 stored=[stored;s1];
end
m1=(dim*(b+m))+b+m;
c1=round(dim^2/2);
yc2=stored(:,c1); % Center Cell
ym2=stored(:,m1); %Mesenchymal
yb2=stored(:,1); %Boundary side

figure
plot(exoc,yc);
hold on
plot(exoc,ym);
hold on
plot(exoc,yb);
ylabel('N-cadherin')
xlabel('TGFB')
legend('Cell Center','Mesenchymal Corner','Boundary','Location','southeast')


%% NCad vs TGFB Case 3
% Testing Case 3 Large epithelial (5E 1M)- With inhibitor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

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

init.mir2=0.0008460221;

init.mask=1;
init.time=[0:1:24*20];% hours


%
n =5; %  dimension  E cells
b=1;  % outer most boundary of no cells
m=1; % middle layer of M cells
dim=n+2*b+2*m; % Dimension of full matrix

exoc=[];
stored=[];
for i=1:81
init.exo=((i-1)/20); %Exogenous concentration of TGF-? 0?5 unit ******Can be changed to 3.0
exoc=[exoc; init.exo];
tic
disp(['Running Cell Simulation with exo: ' num2str(init.exo)])
f1=Multicell(init,n,b,m);
toc
s1=f1(end,end-dim^2+1:end);
 stored=[stored;s1];
end
m1=(dim*(b+m))+b+m;
c1=round(dim^2/2);
yc3=stored(:,c1); % Center Cell
ym3=stored(:,m1); %Mesenchymal
yb3=stored(:,1); %Boundary side

figure
plot(exoc,yc);
hold on
plot(exoc,ym);
hold on
plot(exoc,yb);
ylabel('N-cadherin')
xlabel('TGFB')
legend('Cell Center','Mesenchymal Corner','Boundary','Location','southeast')



%% NCad vs TGFB Case 4
% Testing Case 4 Large epithelial (1E 5M)- With inhibitor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

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

init.mir2=0.0008460221;

init.mask=1;
init.time=[0:1:24*20];% hours


%
n =1; %  dimension  E cells
b=1;  % outer most boundary of no cells
m=5; % middle layer of M cells
dim=n+2*b+2*m; % Dimension of full matrix

exoc=[];
stored=[];
for i=1:81
init.exo=((i-1)/20); %Exogenous concentration of TGF-? 0?5 unit ******Can be changed to 3.0
exoc=[exoc; init.exo];
tic
disp(['Running Cell Simulation with exo: ' num2str(init.exo)])
f1=Multicell(init,n,b,m);
toc
s1=f1(end,end-dim^2+1:end);
 stored=[stored;s1];
end
m1=(dim*(b+m))+b+m;
c1=round(dim^2/2);
yc4=stored(:,c1); % Center Cell
ym4=stored(:,m1); %Mesenchymal
yb4=stored(:,1); %Boundary side

figure
plot(exoc,yc);
hold on
plot(exoc,ym);
hold on
plot(exoc,yb);
ylabel('N-cadherin')
xlabel('TGFB')
legend('Cell Center','Mesenchymal Corner','Boundary','Location','southeast')

%%


%% All cases

subplot(2,3,1)
t2=linspace(0,4,81);
plot(t2,x1);
hold on
plot(t2,x2);
ylabel('N-cadherin')
xlabel('Exogenous TGFB')
legend('Epithelial','Mesenchymal','Location','southeast')
title('Single Cell- No Inhibitor');

subplot(2,3,4)
t2=linspace(0,4,81);
plot(t2,x3);
hold on
plot(t2,x4);
ylabel('N-cadherin')
xlabel('Exogenous TGFB')
legend('Epithelial','Mesenchymal','Location','southeast')
title('Single Cell- Inhibitor');

subplot(2,3,2)
plot(exoc,yc1);
hold on
plot(exoc,ym1);

ylabel('N-cadherin')
xlabel('Exogenous TGFB')
legend('Epithelial Center','Mesenchymal Corner','Location','southeast')
title('Multicell 5E 1M- No Inhibitor');


subplot(2,3,3)
plot(exoc,yc2);
hold on
plot(exoc,ym2);

ylabel('N-cadherin')
xlabel('Exogenous TGFB')
legend('Epithelial Center','Mesenchymal Corner','Location','southeast')
title('Multicell 1E 5M- No Inhibitor');

subplot(2,3,5)
plot(exoc,yc3);
hold on
plot(exoc,ym3);

ylabel('N-cadherin')
xlabel('Exogenous TGFB')
legend('Epithelial Center','Mesenchymal Corner','Location','southeast')
title('Multicell 5E 1M- Inhibitor');


subplot(2,3,6)
plot(exoc,yc4);
hold on
plot(exoc,ym4);
ylabel('N-cadherin')
xlabel('Exogenous TGFB')
legend('Epithelial Center','Mesenchymal Corner','Location','southeast')
title('Multicell 1E 5M- Inhibitor');

%%
load('11sims.mat');

tim=linspace(0,3,31);

plot(tim,Cs1{1,:,1})

