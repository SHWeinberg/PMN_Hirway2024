function make_movie_or_snapshots_v3(name, movie_or_snap, mcs_vec)
% name - filename (without mat)
% movie_or_snap - string 'movie' or 'snap' to set whether to create an
% animation (.mp4) or series of imagesc (.eps)
% mcs_vec - vector that species time points for image creation.  Note that
% for movie making, the vector should be mcs_start:mcs_stop. 
% V3 has phenotype videos for each cell type

%LOAD PARAMETERS
REF=0.1; AMPL=6;
SPACING=2; % show strains every # elements

% load('sim_params.mat');
[cmap,simout]=load_data(name);

if ~isfield(simout{1},'params')
    simout{1}.params = simout{1}.def;
end

def=simout{1}.params;
%     [name,ctag,conn_mat,NRc,jfx,jfy,jxnet,jynet,fx,fy,...
%         fxnet,fynet,pstrain]=parse_data(sim_data{g});
name = simout{1}.name;
if strcmp(name(end-3:end),'.mat')
    name = name(1:end-4);
end
ctag = simout{1}.ctag;
conn_mat = simout{1}.conn_mat;
%     NRc = simout{g}.NRc;
jfx = simout{1}.jfx;
jfy = simout{1}.jfy;
%     jxnet = simout{g}.jxnet;
%     jynet = simout{g}.jynet;
fx = simout{1}.fx;
fy = simout{1}.fy;
%     fxnet = simout{g}.fxnet;
%     fynet = simout{g}.fynet;
pstrain = simout{1}.pstrain;
pstrainx = simout{1}.pstrainx;
pstrainy = simout{1}.pstrainy;
phenotype = simout{1}.phenotype;
celltype = simout{1}.pop_index;

%INITIALIZE GLOBALS
NVX=def.NVX; NVY=def.NVY;
NVX2=SPACING*NVX; NVY2=SPACING*NVY;
NNX=def.NNX; NNY=def.NNY; NN=def.NN;

NRINC = length(mcs_vec);


% calculate max/min strain magnitudes
strsc_all = pstrain(:,mcs_vec)*1e-6;
strx_all = strsc_all.*pstrainx(:,mcs_vec)/1e3;
stry_all = strsc_all.*pstrainy(:,mcs_vec)/1e3;
strmag_all = sqrt(strx_all.^2 + stry_all.^2);
str_range = [min(strmag_all(:)) max(strmag_all(:))];

%INITIALIZE FIGURES
[cellim,cellax]=cellmap_init();
[cfig2,cellax2]=cellmap_init2();
[cfig3,cellax2]=cellmap_init2();
[cfig4,cellax2]=cellmap_init2();
[cfig5,cellax2]=cellmap_init2();
[cfig6,cellax2]=cellmap_init2();

[s,strfig]=strainmap_init();
[cntrax,cntr,q,allfig,allfig2]=contourmap_init();
% [networkax,networkfig]=networkgraph_init();
%     [b,barfig]=forcebars_init(NRc{1});

%GENERATE FIGURES
Mc=cell(NRINC,1);
Ms=cell(NRINC,1);
Mg=cell(NRINC,1);
Mb=cell(NRINC,1);
Ma=cell(NRINC,1);
Ma2=cell(NRINC,1);
Mc2=cell(NRINC,1);
Mt=cell(NRINC,1);
MK=cell(NRINC,1);
MS=cell(NRINC,1);
MT=cell(NRINC,1);
count = 0;
for c=mcs_vec
    count = count+1;
    sprintf('ITERATION %f of %f',c,NRINC)
    strfield = map_strain(s,cntr,cmap.depth,pstrain(:,c), pstrainx(:,c), pstrainy(:,c));
    im = map_cells(cellim, cmap,ctag(:,c),jfx(:,c),jfy(:,c));
    map_cells2(cfig2, cmap,ctag(:,c),phenotype{c}, im);
    map_cells3(cfig3, cmap,ctag(:,c),celltype{c}, im);
     map_cellsK(cfig4, cmap,ctag(:,c),phenotype{c},celltype{c}, im);
     map_cellsS(cfig5, cmap,ctag(:,c),phenotype{c},celltype{c}, im);
     map_cellsT(cfig6, cmap,ctag(:,c),phenotype{c},celltype{c}, im);
    q = map_forces(cntrax,q,fx(:,c),fy(:,c));
   % networkgraph(networkax,conn_mat{c}, networkfig);
    cell_strain_traction(allfig, allfig2, strfield, im, q, c, str_range);
    %         forcebar(b,jxnet{c},jynet{c},fxnet{c},fynet{c},NRc{c});
    
     Mcols={'CellTags', ...
            'StrainMap', ...
            'NetworkGraph', ...
            'ContourMap', ...
            'Cells_Forces', ...
            'Cells_Forces_norm', ...
            'Cell_Phenotype',...
            'Cell_Types',...
            'KType',...
            'SType',...
            'TType'};
        switch movie_or_snap
            case 'snap'
                saveas(cellax, [Mcols{1},'_',name,'_',num2str(c)],'epsc');
                saveas(strfig, [Mcols{2},'_',name,'_',num2str(c)],'epsc');
                saveas(networkfig, [Mcols{3},'_',name,'_',num2str(c)],'epsc');
                saveas(cntrax, [Mcols{4},'_',name,'_',num2str(c)],'epsc');
                saveas(allfig, [Mcols{5},'_',name,'_',num2str(c)],'epsc');
                saveas(allfig2, [Mcols{6},'_',name,'_',num2str(c)],'epsc');
                saveas(cfig2, [Mcols{7},'_',name,'_',num2str(c)],'epsc');
                saveas(cfig3, [Mcols{8},'_',name,'_',num2str(c)],'epsc');
                saveas(cfig4, [Mcols{9},'_',name,'_',num2str(c)],'epsc');
                saveas(cfig5, [Mcols{10},'_',name,'_',num2str(c)],'epsc');
                saveas(cfig6, [Mcols{11},'_',name,'_',num2str(c)],'epsc');
            case 'movie'
                Mc{count}=getframe(cellax);
              %  Ms{count}=getframe(strfig);
%                 Mg{count}=getframe(networkfig);
%                 Mb{count}=getframe(cntrax);
%                 Ma{count}=getframe(allfig);
                Ma2{count}=getframe(allfig2);
                Mc2{count}=getframe(cfig2);
                Mt{count}=getframe(cfig3);
                MK{count}=getframe(cfig4);
                MS{count}=getframe(cfig5);
                MT{count}=getframe(cfig6);
                
        end
end
if strcmp(movie_or_snap,'movie')
%     Mdata=[{Mc}, {Ms}, {Mg}, {Mb}, {Ma},{Ma2}, {Mc2}];

    %SAVE MOVIE
    write_mov(name,Mc,Mcols{1});
%     write_mov(name,Ms,Mcols{2});
%     write_mov(name,Mg,Mcols{3});
%     write_mov(name,Mb,Mcols{4});
%     write_mov(name,Ma,Mcols{5});
    write_mov(name,Ma2,Mcols{6});
    write_mov(name,Mc2,Mcols{7});
    write_mov(name,Mt,Mcols{8});
    write_mov(name,MK,Mcols{9});
    write_mov(name,MS,Mcols{10});
    write_mov(name,MT,Mcols{11});


end



%_________________________________________________________________________%
    function [cmap,simout]=load_data(name)
        %         listing=dir('sim_data_*');
        %         simout=cell(length(listing),1);
        %         sim_data=cell(length(listing),1);
        %         for ln=1:length(listing)
        simout{1}=load(name);
%         
%         if strfind(name,'atag')
%             junc_name = [name,'_','junction'];
%             juncout = load(junc_name);
%             
%             jfx = zeros(simout{1}.params.NV, simout{1}.params.NRINC);
%             jfy = zeros(simout{1}.params.NV, simout{1}.params.NRINC);
%             
%             for i = 1:simout{1}.params.NRINC
%                 xvox_i = juncout.xjunc_vox{i};
%                 yvox_i = juncout.yjunc_vox{i};
%                 try
%                     ind = sub2ind([simout{1}.params.NVY simout{1}.params.NVX], ...
%                         yvox_i, xvox_i);
%                 catch
%                     1;
%                 end
%                 jfx(ind,i) = juncout.jx_vox_norm{i};
%                 jfy(ind,i) = juncout.jy_vox_norm{i};
%                 
%             end
%             simout{1}.jfx = jfx;
%             simout{1}.jfy = jfy;
%         end
        
        %             sim_data{ln}=simout{ln}.sim_data;
        %         end
        
        % black=129 gray=130 green=131 white=132 blue=133
        colorcols={
            'depth', ...
            'black', ...
            'gray', ...
            'green', ...
            'white', ...
            'blue'};
        cmap.depth=256; % color map bit depth
        cmapdata=cell(1,6);
        for k=0:5
            cmapdata{k+1}=cmap.depth+k;
        end
        cmap=cell2struct(cmapdata,colorcols,2);
    end

    function [name,ctag,conn_mat,NRc,jfx,jfy,jxnet,jynet,fx,fy,fxnet,fynet,pstrain]=parse_data(sim_data)
        name=sim_data.Name;
        ctag={sim_data.Cells.CellLocation};
        conn_mat={sim_data.Cells.ConnectivityMatrix};
        NRc={sim_data.Cells.CellCount};
        
        jf=sim_data.FE.Junction.Cell;
        jfx={jf.XData}; jfy={jf.YData};
        
        jfnet=sim_data.FE.Junction.Net;
        jxnet={jfnet.XData}; jynet={jfnet.YData};
        
        tf=sim_data.FE.Traction.Cell;
        fx={tf.XData}; fy={tf.YData};
        
        tfnet=sim_data.FE.Traction.Net;
        fxnet={tfnet.XData}; fynet={tfnet.YData};
        
        pstrain=sim_data.FE.FEStrain;
    end

%_________________________________________________________________________%

    function [cellim,cellax]=cellmap_init
        cfig=figure('name','Cell Tags');
        set(cfig,'position', [100 200+2*NVY2 3*NVX2 3*NVY2])
        
        bigfield = 260*ones(NVX2,NVY2);
        cellax=axes(cfig);
        cellim=image(cellax,bigfield);
        colormap(cellax,parula2)
        axis(cellax,'image'), axis(cellax,'xy'), axis(cellax,'off')
        
        title(cellax,'Cellular Potts Visualization')
    end


    function [cfig,cellax]=cellmap_init2
        cfig=figure('name','Cell Phenotype');
        set(cfig,'position', [100 200+2*NVY2 3*NVX2 3*NVY2])
        
        bigfield = 260*zeros(NVX2,NVY2);
        cellax=axes(cfig);
        cellim2=image(cellax,bigfield);
        colormap(cellax,jet); caxis([0 1]);
        axis(cellax,'image'), axis(cellax,'xy'), axis(cellax,'off')
        
        title(cellax,'Cellular Potts Phenotype Visualization')
    end

function [cfig,cellax]=cellmap_init3
        cfig=figure('name','Cell Types');
        set(cfig,'position', [100 200+2*NVY2 3*NVX2 3*NVY2])
        
        bigfield = 260*zeros(NVX2,NVY2);
        cellax=axes(cfig);
        cellim2=image(cellax,bigfield);
        colormap(cellax,jet); caxis([0 1]);
        axis(cellax,'image'), axis(cellax,'xy'), axis(cellax,'off')
        
        title(cellax,'Cellular Potts Cell Types Visualization')
    end

    function [s,strfig]=strainmap_init
        strfig=figure('name','Strain Map');
        strax=axes(strfig);
        set(strfig,'position',[100+1.1*2*NVX2 200+2*NVY2 3*NVX2 2.5*NVY2])
        strfield=zeros(NVX2,NVY2);
        X=1:NVX2; Y=1:NVY2;
        s=surf(X,Y,strfield);
        s.FaceColor='interp';
        s.FaceAlpha=0.8;
        s.EdgeAlpha=0.7;
        
        xyrot=get(strax,'view');
        xlabel(strax,'X (\mum)','FontSize',12,'FontWeight','bold','rotation',xyrot(2))
        ylabel(strax,'Y (\mum)','FontSize',12,'FontWeight','bold','rotation',xyrot(1))
        zlabel(strax,'\epsilon','FontSize',14,'FontWeight','bold')
        
        strax.XTickMode='manual';
        strax.YTickMode='manual';
        surfticks=[0 50 100 150 200];
        surfticklbls={'0' '125' '250' '375' '500'};
        strax.XTick=surfticks;
        strax.YTick=surfticks;
        strax.XTickLabel=surfticklbls;
        strax.YTickLabel=surfticklbls;
        
        colormap(strax,'parula'), colorbar(strax)
        title('Finite Element Strain')
    end

    function [cntrax,cntr,q, allfig,allfig2]=contourmap_init()
        cntrfig=figure('name','Strain Map');
        cntrax=axes(cntrfig);
        set(cntrfig,'position',[100+1.1*2*NVX2 200+2*NVY2 3*NVX2 2.5*NVY2])
        axis(cntrax,'image')
        cntrax.XTickMode='manual';
        cntrax.YTickMode='manual';
        cntrax.XTick=[];
        cntrax.YTick=[];
        
        
        colormap(cntrax,'parula')
        colorbar(cntrax)
        
        cntrfield=zeros(NVX2,NVY2);
        X=1:NVX2; Y=1:NVY2;
        [~,cntr]=contour(X,Y,cntrfield);
        cntr.LineWidth=2;
        cntr.Fill='on';
        
        hold(cntrax,'on')
        q=quiver(cntrax,[],[],[],[],'w');
        q.LineWidth=1;
        q.AutoScaleFactor=1.5;
        hold(cntrax,'off')
        title('FE strain with traction force overlay')
        
        allfig=figure('name','Cell positions and forces');
        set(allfig,'position',[50+1.1*2*NVX2 200+2*NVY2 3*NVX2 2.5*NVY2],'color','w')
        
        allfig2=figure('name','Cell positions and forces - normalized');
        set(allfig2,'position',[60+1.1*2*NVX2 200+2*NVY2 3*NVX2 2.5*NVY2],'color','w')
        
    end

    function [networkax,networkfig]=networkgraph_init()
        networkfig=figure('name','Network Graph');
        set(networkfig,'position',[200 120 3*NVX2 2.5*NVY2])
        networkax=axes(networkfig);
        colormap(networkax,'parula')
        
        networkax.XTickMode='manual';
        networkax.YTickMode='manual';
        networkax.XTick=[];
        networkax.YTick=[];
        networkax.Box='on';
        axis(networkax,'image');
        
        title(networkax,'Cellular Network Connectivity')
        drawnow;
    end

    function [b,barfig]=forcebars_init(NRc)
        barfig=figure('name','Force Comparison');
        set(barfig,'position',[200+1.1*2*NVX2 120 3*NVX2 3*NVY2])
        barax=axes(barfig);
        
        X=1:NRc; Y=[zeros(NRc, 1), zeros(NRc, 1)];
        b=bar(barax,X,Y);
        
        title(barax,'FMA net junction and traction forces per cell')
        xlabel(barax,'Cell Number','FontWeight','Bold')
        ylabel(barax,'Force Magnitude (nN)','FontWeight','Bold')
        
        lg=legend(barax,[{'Junction Force'}, {'Traction Force'}]);
        lg.FontWeight='bold';
        lg.Location='best';
    end

%_________________________________________________________________________%

%STRAIN DISPLAY
    function strfield = map_strain(s,cntr,m,pstrain, pstrainx, pstrainy)
        %SUBSTRATE STRAIN
        strfield=zeros(NVX2,NVY2);
        strsc = pstrain(:,1)*1E-6;
        %         strsc = AMPL*str/REF;
        strx=zeros(NVX,NVY); stry=zeros(NVX,NVY);
        netstr = @(x,y) (sqrt(x^2 + y^2));
        for ey=1:NVY
            for ex=1:NVX
                e = ex+(ey-1)*NVX;
                strx(ex,ey)=strsc(e)*pstrainx(e)/1000;
                stry(ex,ey)=strsc(e)*pstrainy(e)/1000;
                strfield(2*ex-1,2*ey-1)=netstr(strx(ex,ey), stry(ex,ey));
                strfield(2*ex-1,2*ey  )=netstr(strx(ex,ey), stry(ex,ey));
                strfield(2*ex  ,2*ey-1)=netstr(strx(ex,ey), stry(ex,ey));
                strfield(2*ex  ,2*ey  )=netstr(strx(ex,ey), stry(ex,ey));
            end
        end
        cmin=min(strfield(:)); cmax=max(strfield(:));
        bigfield = min(m,(m-1)*(strfield-cmin)/(cmax-cmin)+1);
        set(s,'zdata',bigfield);
        set(cntr,'zdata',bigfield);
    end

%CELL DISPLAY
    function im = map_cells(cellim,cmap,ctag,jfx,jfy)
        %CELL FIGURE
        bigctags = zeros(NVX2,NVY2);
        bigjmag = zeros(NVX2,NVY2);
        jmag=sqrt(jfx.^2+jfy.^2);
        jmagmin=min(jmag); jmagmax=max(jmag);
        jmagsc=(cmap.depth-1)*(jmag-jmagmin)/(jmagmax-jmagmin)+1;
        jmagsc(jmag==jmagmin)=0;
        for ey = 1:NVY
            for ex = 1:NVX
                e = ex+(ey-1)*NVX;
                bigctags(ey*2-1,ex*2-1)=ctag(e);
                bigctags(ey*2-1,ex*2  )=ctag(e);
                bigctags(ey*2  ,ex*2-1)=ctag(e);
                bigctags(ey*2  ,ex*2  )=ctag(e);
                
                bigjmag(ey*2-1,ex*2-1)=jmagsc(e);
                bigjmag(ey*2-1,ex*2  )=jmagsc(e);
                bigjmag(ey*2  ,ex*2-1)=jmagsc(e);
                bigjmag(ey*2  ,ex*2  )=jmagsc(e);
            end
        end
        %         bigfield=get(cellim,'CData');
        bigfield=260*ones(NVX2,NVY2);
        for ey=2:NVY2
            for ex=2:NVX2
                if (bigctags(ey,ex)) 
                   bigfield(ey,ex) = 259;  % 258 = red, 259 = gray
                    
                end
                if (bigctags(ey,ex) ~= bigctags(ey-1,ex)) || (bigctags(ey,ex) ~= bigctags(ey,ex-1))
                    if bigjmag(ey,ex)~=0
                        bigfield(ey,ex)=bigjmag(ey,ex);
                    else
                        bigfield(ey,ex) = 257;  
                        
                    end
                end
            end
        end
        im = zeros(size(bigfield)); 
       im(bigfield<259) = 1; 
       
        set(cellim,'CData',bigfield)
        %         set(cntrax,'cdata',bigfield);
        drawnow
    end

 function map_cells2(cfig,cmap,ctag,pmag, im)
        figure(cfig); clf
        
        ax1 = axes('parent',cfig); ax2 = axes('parent',cfig);
     %CELL FIGURE
        bigctags = zeros(NVX2,NVY2);
        bigpmag = zeros(NVX2,NVY2);

        for ey = 1:NVY
            for ex = 1:NVX
                e = ex+(ey-1)*NVX;
                bigctags(ey*2-1,ex*2-1)=ctag(e);
                bigctags(ey*2-1,ex*2  )=ctag(e);
                bigctags(ey*2  ,ex*2-1)=ctag(e);
                bigctags(ey*2  ,ex*2  )=ctag(e);

                if ctag(e)
                    bigpmag(ey*2-1,ex*2-1)=pmag(ctag(e));
                    bigpmag(ey*2-1,ex*2  )=pmag(ctag(e));
                    bigpmag(ey*2  ,ex*2-1)=pmag(ctag(e));
                    bigpmag(ey*2  ,ex*2  )=pmag(ctag(e));
                end
            end
        end
        ap = ones(size(bigpmag));
        ap(bigctags==0) = 0;
        ap(im==1) = 0;
        % COULD CHANGE BORDER COLOR HERE- SUH 10/29/21
       %ap(im==1)= 0.25;
        axes(ax2);
        imagesc(bigpmag); alpha(ap);
        caxis([0 1]); colormap(jet); hold on;
        I = imshow(1-im,'parent',ax1);   
        set(I,'alphadata',im);
       
          title(['MCS ', num2str(c)]);
        set(ax1,'fontsize',26,'ydir','normal','xlim',[.5 NVX2+.5],'ylim',[.5 NVY2+.5]);
        set(ax2,'fontsize',26,'ydir','normal','xlim',[.5 NVX2+.5],'ylim',[.5 NVY2+.5]);
        axis(ax1,'image'); 
        axis(ax2,'image');
        axis(ax1,'xy'); axis(ax2,'xy');
        axis(ax1,'off'); axis(ax2,'off');
        
        %         bigfield=get(cellim,'CData');
%         bigfield=260*ones(NVX2,NVY2);
%         for ey=2:NVY2
%             for ex=2:NVX2
%                 if (bigctags(ey,ex))  %
%                     bigfield(ey,ex) = bigpmag(ey,ex); %259;  % 258 = red, 259 = gray
%                 end
%                 if (bigctags(ey,ex) ~= bigctags(ey-1,ex)) || (bigctags(ey,ex) ~= bigctags(ey,ex-1))
%                     %if bigjmag(ey,ex)~=0
%                     %    bigfield(ey,ex)= 257; bigjmag(ey,ex);
%                     %else
%                         bigfield(ey,ex) = 257;
%                     %end
%                 end
%             end
%         end
%         im = zeros(size(bigfield)); im(bigfield<259) = 1;
%         set(cellim,'CData',bigfield)
        %         set(cntrax,'cdata',bigfield);
        drawnow
 end


function map_cells3(cfig,cmap,ctag,pmag, im) % For cell types***************** SUH 11/17/21
        figure(cfig); clf
        
        ax1 = axes('parent',cfig); ax2 = axes('parent',cfig);
     %CELL FIGURE
        bigctags = zeros(NVX2,NVY2);
        bigpmag = zeros(NVX2,NVY2);

        for ey = 1:NVY
            for ex = 1:NVX
                e = ex+(ey-1)*NVX;
                bigctags(ey*2-1,ex*2-1)=ctag(e);
                bigctags(ey*2-1,ex*2  )=ctag(e);
                bigctags(ey*2  ,ex*2-1)=ctag(e);
                bigctags(ey*2  ,ex*2  )=ctag(e);

                if ctag(e)
                    bigpmag(ey*2-1,ex*2-1)=pmag(ctag(e));
                    bigpmag(ey*2-1,ex*2  )=pmag(ctag(e));
                    bigpmag(ey*2  ,ex*2-1)=pmag(ctag(e));
                    bigpmag(ey*2  ,ex*2  )=pmag(ctag(e));
                end
            end
        end
        ap = ones(size(bigpmag));
        ap(bigctags==0) = 0;
        ap(im==1) = 0;
        % COULD CHANGE BORDER COLOR HERE- SUH 10/29/21
       %ap(im==1)= 0.25;
        axes(ax2);
        imagesc(bigpmag); alpha(ap);
        caxis([1 3]); colormap(jet); hold on;
        I = imshow(1-im,'parent',ax1);   
        set(I,'alphadata',im);
       
          title(['MCS ', num2str(c)]);
        set(ax1,'fontsize',26,'ydir','normal','xlim',[.5 NVX2+.5],'ylim',[.5 NVY2+.5]);
        set(ax2,'fontsize',26,'ydir','normal','xlim',[.5 NVX2+.5],'ylim',[.5 NVY2+.5]);
        axis(ax1,'image'); 
        axis(ax2,'image');
        axis(ax1,'xy'); axis(ax2,'xy');
        axis(ax1,'off'); axis(ax2,'off');
        
        %         bigfield=get(cellim,'CData');
%         bigfield=260*ones(NVX2,NVY2);
%         for ey=2:NVY2
%             for ex=2:NVX2
%                 if (bigctags(ey,ex))  %
%                     bigfield(ey,ex) = bigpmag(ey,ex); %259;  % 258 = red, 259 = gray
%                 end
%                 if (bigctags(ey,ex) ~= bigctags(ey-1,ex)) || (bigctags(ey,ex) ~= bigctags(ey,ex-1))
%                     %if bigjmag(ey,ex)~=0
%                     %    bigfield(ey,ex)= 257; bigjmag(ey,ex);
%                     %else
%                         bigfield(ey,ex) = 257;
%                     %end
%                 end
%             end
%         end
%         im = zeros(size(bigfield)); im(bigfield<259) = 1;
%         set(cellim,'CData',bigfield)
        %         set(cntrax,'cdata',bigfield);
        drawnow
end

function map_cellsK(cfig,cmap,ctag,pmag,tmag, im)
        figure(cfig); clf
        Ktype=zeros(size(pmag));
        Ktype(find(tmag==1))=pmag(find(tmag==1));
        Kindex=find(tmag~=1);
        NonK=ismember(ctag,Kindex);
        ctag(NonK)=0;
         ax1 = axes('parent',cfig); ax2 = axes('parent',cfig);
     %CELL FIGURE- K
    
        bigctags = zeros(NVX2,NVY2);
        bigpmag = zeros(NVX2,NVY2);

        for ey = 1:NVY
            for ex = 1:NVX
                e = ex+(ey-1)*NVX;
                bigctags(ey*2-1,ex*2-1)=ctag(e);
                bigctags(ey*2-1,ex*2  )=ctag(e);
                bigctags(ey*2  ,ex*2-1)=ctag(e);
                bigctags(ey*2  ,ex*2  )=ctag(e);

                if ctag(e)
                    bigpmag(ey*2-1,ex*2-1)=Ktype(ctag(e));
                    bigpmag(ey*2-1,ex*2  )=Ktype(ctag(e));
                    bigpmag(ey*2  ,ex*2-1)=Ktype(ctag(e));
                    bigpmag(ey*2  ,ex*2  )=Ktype(ctag(e));
                end
            end
        end
        ap = ones(size(bigpmag));
        ap(bigctags==0) = 0;
        ap(im==1) = 0;
        % COULD CHANGE BORDER COLOR HERE- SUH 10/29/21
       %ap(im==1)= 0.25;
       
        axes(ax2); %commented out
        imagesc(bigpmag); alpha(ap); %image being set
    
        
        caxis([0 1]); colormap(jet); hold on;
        I = imshow(1-im,'parent',ax1);   
        set(I,'alphadata',im);
       
          title(['MCS ', num2str(c)]);
        set(ax1,'fontsize',26,'ydir','normal','xlim',[.5 NVX2+.5],'ylim',[.5 NVY2+.5]);
        set(ax2,'fontsize',26,'ydir','normal','xlim',[.5 NVX2+.5],'ylim',[.5 NVY2+.5]);
        axis(ax1,'image'); 
        axis(ax2,'image');
        axis(ax1,'xy'); axis(ax2,'xy');
        axis(ax1,'off'); axis(ax2,'off');
        
        %         bigfield=get(cellim,'CData');
%         bigfield=260*ones(NVX2,NVY2);
%         for ey=2:NVY2
%             for ex=2:NVX2
%                 if (bigctags(ey,ex))  %
%                     bigfield(ey,ex) = bigpmag(ey,ex); %259;  % 258 = red, 259 = gray
%                 end
%                 if (bigctags(ey,ex) ~= bigctags(ey-1,ex)) || (bigctags(ey,ex) ~= bigctags(ey,ex-1))
%                     %if bigjmag(ey,ex)~=0
%                     %    bigfield(ey,ex)= 257; bigjmag(ey,ex);
%                     %else
%                         bigfield(ey,ex) = 257;
%                     %end
%                 end
%             end
%         end
%         im = zeros(size(bigfield)); im(bigfield<259) = 1;
%         set(cellim,'CData',bigfield)
        %         set(cntrax,'cdata',bigfield);
        drawnow
end

function map_cellsS(cfig,cmap,ctag,pmag,tmag, im)
        figure(cfig); clf
      
        Stype=zeros(size(pmag));
         Stype(find(tmag==2))=pmag(find(tmag==2));
          Sindex=find(tmag~=2);
        NonS=ismember(ctag,Sindex);
        ctag(NonS)=0;
        
         ax1 = axes('parent',cfig); ax2 = axes('parent',cfig);
     %CELL FIGURE- S
    
        bigctags = zeros(NVX2,NVY2);
        bigpmag = zeros(NVX2,NVY2);

        for ey = 1:NVY
            for ex = 1:NVX
                e = ex+(ey-1)*NVX;
                bigctags(ey*2-1,ex*2-1)=ctag(e);
                bigctags(ey*2-1,ex*2  )=ctag(e);
                bigctags(ey*2  ,ex*2-1)=ctag(e);
                bigctags(ey*2  ,ex*2  )=ctag(e);

                if ctag(e)
                    bigpmag(ey*2-1,ex*2-1)=Stype(ctag(e));
                    bigpmag(ey*2-1,ex*2  )=Stype(ctag(e));
                    bigpmag(ey*2  ,ex*2-1)=Stype(ctag(e));
                    bigpmag(ey*2  ,ex*2  )=Stype(ctag(e));
                end
            end
        end
        ap = ones(size(bigpmag));
        ap(bigctags==0) = 0;
        ap(im==1) = 0;
        % COULD CHANGE BORDER COLOR HERE- SUH 10/29/21
       %ap(im==1)= 0.25;
       
        axes(ax2); %commented out
        imagesc(bigpmag); alpha(ap); %image being set
    
        
        caxis([0 1]); colormap(jet); hold on;
        I = imshow(1-im,'parent',ax1);   
        set(I,'alphadata',im);
       
          title(['MCS ', num2str(c)]);
        set(ax1,'fontsize',26,'ydir','normal','xlim',[.5 NVX2+.5],'ylim',[.5 NVY2+.5]);
        set(ax2,'fontsize',26,'ydir','normal','xlim',[.5 NVX2+.5],'ylim',[.5 NVY2+.5]);
        axis(ax1,'image'); 
        axis(ax2,'image');
        axis(ax1,'xy'); axis(ax2,'xy');
        axis(ax1,'off'); axis(ax2,'off');
        
        %         bigfield=get(cellim,'CData');
%         bigfield=260*ones(NVX2,NVY2);
%         for ey=2:NVY2
%             for ex=2:NVX2
%                 if (bigctags(ey,ex))  %
%                     bigfield(ey,ex) = bigpmag(ey,ex); %259;  % 258 = red, 259 = gray
%                 end
%                 if (bigctags(ey,ex) ~= bigctags(ey-1,ex)) || (bigctags(ey,ex) ~= bigctags(ey,ex-1))
%                     %if bigjmag(ey,ex)~=0
%                     %    bigfield(ey,ex)= 257; bigjmag(ey,ex);
%                     %else
%                         bigfield(ey,ex) = 257;
%                     %end
%                 end
%             end
%         end
%         im = zeros(size(bigfield)); im(bigfield<259) = 1;
%         set(cellim,'CData',bigfield)
        %         set(cntrax,'cdata',bigfield);
        drawnow
end

function map_cellsT(cfig,cmap,ctag,pmag,tmag, im)
        figure(cfig); clf
     
        Ttype=zeros(size(pmag));
 
          Ttype(find(tmag==3))=pmag(find(tmag==3));
         Tindex=find(tmag~=3);
        NonT=ismember(ctag,Tindex);
        ctag(NonT)=0;
          
         ax1 = axes('parent',cfig); ax2 = axes('parent',cfig);
     %CELL FIGURE- T
    
        bigctags = zeros(NVX2,NVY2);
        bigpmag = zeros(NVX2,NVY2);

        for ey = 1:NVY
            for ex = 1:NVX
                e = ex+(ey-1)*NVX;
                bigctags(ey*2-1,ex*2-1)=ctag(e);
                bigctags(ey*2-1,ex*2  )=ctag(e);
                bigctags(ey*2  ,ex*2-1)=ctag(e);
                bigctags(ey*2  ,ex*2  )=ctag(e);

                if ctag(e)
                    bigpmag(ey*2-1,ex*2-1)=Ttype(ctag(e));
                    bigpmag(ey*2-1,ex*2  )=Ttype(ctag(e));
                    bigpmag(ey*2  ,ex*2-1)=Ttype(ctag(e));
                    bigpmag(ey*2  ,ex*2  )=Ttype(ctag(e));
                end
            end
        end
        ap = ones(size(bigpmag));
        ap(bigctags==0) = 0;
        ap(im==1) = 0;
        % COULD CHANGE BORDER COLOR HERE- SUH 10/29/21
       %ap(im==1)= 0.25;
       
        axes(ax2); %commented out
        imagesc(bigpmag); alpha(ap); %image being set
    
        
        caxis([0 1]); colormap(jet); hold on;
        I = imshow(1-im,'parent',ax1);   
        set(I,'alphadata',im);
       
          title(['MCS ', num2str(c)]);
        set(ax1,'fontsize',26,'ydir','normal','xlim',[.5 NVX2+.5],'ylim',[.5 NVY2+.5]);
        set(ax2,'fontsize',26,'ydir','normal','xlim',[.5 NVX2+.5],'ylim',[.5 NVY2+.5]);
        axis(ax1,'image'); 
        axis(ax2,'image');
        axis(ax1,'xy'); axis(ax2,'xy');
        axis(ax1,'off'); axis(ax2,'off');
        
        %         bigfield=get(cellim,'CData');
%         bigfield=260*ones(NVX2,NVY2);
%         for ey=2:NVY2
%             for ex=2:NVX2
%                 if (bigctags(ey,ex))  %
%                     bigfield(ey,ex) = bigpmag(ey,ex); %259;  % 258 = red, 259 = gray
%                 end
%                 if (bigctags(ey,ex) ~= bigctags(ey-1,ex)) || (bigctags(ey,ex) ~= bigctags(ey,ex-1))
%                     %if bigjmag(ey,ex)~=0
%                     %    bigfield(ey,ex)= 257; bigjmag(ey,ex);
%                     %else
%                         bigfield(ey,ex) = 257;
%                     %end
%                 end
%             end
%         end
%         im = zeros(size(bigfield)); im(bigfield<259) = 1;
%         set(cellim,'CData',bigfield)
        %         set(cntrax,'cdata',bigfield);
        drawnow
 end


%FORCE VECTORS
    function q = map_forces(cntrax,q,fx,fy)
        hold(cntrax,'on')
        fxsc = AMPL*fx/REF; fysc = AMPL*fy/REF;
        xpos=zeros(NN,1); ypos=zeros(NN,1);
        forcex=zeros(NNX,1); forcey=zeros(NNY,1);
        for ey=1:SPACING:NNY
            for ex=1:SPACING:NNX
                e = ex+(ey-1)*NNX;
                xpos(e)=ex*2;
                ypos(e)=ey*2;
                forcex(e) = fxsc(e);
                forcey(e) = fysc(e);
            end
        end
        set(q,'xdata',xpos,'ydata',ypos,'udata',forcex,'vdata',forcey)
        hold(cntrax,'off')
        drawnow
    end

%NETWORK GRAPH OF CELL CONNECTIVITY
    function networkgraph(networkax,conn_mat,netfig)
        figure(netfig);
        % CONNECTIVITY GRAPH
        G=graph(conn_mat);
        deg = degree(G);
        G=rmnode(G,find(deg==0));
        deg=deg(find(deg));
        nsizes = 2*sqrt(deg-min(deg)+2);
        
        p=plot(G);
        p.EdgeAlpha=.4;
        p.LineWidth=1.5;
        p.EdgeColor='k';
        p.NodeCData=deg;
        p.MarkerSize=nsizes;
        colorbar(networkax)
        %         layout(p,'force3')
        
        networkax.XTick=[];
        networkax.YTick=[];
        networkax.ZTick=[];
        networkax.Box='on';
        
        title(networkax,'Cellular Network Connectivity')
        drawnow;
    end

%NET FORCE COMPARISON PER CELL
    function forcebar(b,jxnet,jynet,fxnet,fynet,NRc)
        jnet=sqrt(jxnet.^2+jynet.^2);
        fnet=sqrt(fxnet.^2+fynet.^2);
        
        X=1:NRc;
        set(b(1),'xdata',X,'ydata',jnet)
        set(b(2),'xdata',X,'ydata',fnet)
        drawnow
    end

    function write_mov(name,Mdata, pre)
       
       
        vid_file=[pre, '_', name];
        v = VideoWriter(vid_file,'MPEG-4');
        v.Quality=100; v.FrameRate=15;
        open(v)
        for n=1:NRINC
            writeVideo(v,Mdata{n})
        end
        close(v)
        
    end

    function cell_strain_traction(fighandle, fighandle2, strfield, im, q, c, str_range)
        
        X=1:NVX2; Y=1:NVY2;
        
        figure(fighandle); clf
        
        ax1 = axes('parent',fighandle); ax2 = axes('parent',fighandle);
        contour(X,Y,strfield,'parent',ax1,'linewidth',2,'fill','on');
        I = imshow(1-im,'parent',ax2);
        set(I,'alphadata',im);
        colormap(jet);
        hold on;
        qq = quiver(q.XData,q.YData,q.UData,q.VData,'w');
        qq.LineWidth=1;
        qq.AutoScaleFactor=1.5;
        hold off;
        title(['MCS ', num2str(c)]);
        set(ax1,'fontsize',26,'ydir','normal','xlim',[.5 NVX2+.5],'ylim',[.5 NVY2+.5]);
        set(ax2,'fontsize',26,'ydir','normal','xlim',[.5 NVX2+.5],'ylim',[.5 NVY2+.5]);
        axis(ax1,'image'); axis(ax2,'image');
        axis(ax1,'xy'); axis(ax2,'xy');
        axis(ax1,'off'); axis(ax2,'off');
        
        figure(fighandle2); clf
        
        ax1 = axes('parent',fighandle2); ax2 = axes('parent',fighandle2);
        contour(X,Y,strfield,'parent',ax1,'linewidth',2,'fill','on');
        set(ax1,'clim',str_range);
        I = imshow(1-im,'parent',ax2);
        set(I,'alphadata',im);
        colormap(jet);
        hold on;
        qq = quiver(q.XData,q.YData,q.UData,q.VData,'w');
        qq.LineWidth=1;
        qq.AutoScaleFactor=1.5;
        hold off;
        title(['MCS ', num2str(c)]);
        set(ax1,'fontsize',26,'ydir','normal','xlim',[.5 NVX2+.5],'ylim',[.5 NVY2+.5]);
        set(ax2,'fontsize',26,'ydir','normal','xlim',[.5 NVX2+.5],'ylim',[.5 NVY2+.5]);
        axis(ax1,'image'); axis(ax2,'image');
        axis(ax1,'xy'); axis(ax2,'xy');
        axis(ax1,'off'); axis(ax2,'off');
        
    end

end
