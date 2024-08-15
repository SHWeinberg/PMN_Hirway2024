function [ctag,NRc,csize, statevars,parammatrix,phenoparams,cellmod]= ...
    cell_proliferation_v3(ctag,NRc,csize,def, statevars, pdivs, targetvols,parammatrix,cellmod,phenoparams)
NVX=def.NVX;
NVY=def.NVY;
% divthresh=def.divthresh;
%split=cellmod.split;
Amin=2/3;
% Aideal=def.TARGETVOLUME;                        %SS cell area

ncells=NRc;
 phenoparams=phenoparams;
[~,nstate] = size(statevars);

for c=1:ncells                                  %prob for each cell
    cellA=csize(c);
    Aideal = targetvols(c);
    Athresh=Amin*Aideal;                            %minimum area for division
%  PROLIFERATION RATE FOR CELLS BEFORE CELL 93 (so, 1-92) is set to 0   -------------------- ADDED HERE 8/28/21
% if(c<split)
%     pdivs(c)=0;
% end
    if cellA>Athresh
        divthresh = pdivs(c);
        if divthresh>rand
            cellind=find(ctag==c);
            vx=rem(cellind-1,NVX)+1; vy=(cellind-vx)/NVY + 1;
            Cx=mean(vx); Cy=mean(vy);
            Ixx=sum((vy-Cy).^2);
            Iyy=sum((vx-Cx).^2);
            Ixy=-sum((vx-Cx).*(vy-Cy));
            
            lam=0.5*(Ixx+Iyy)+0.5*sqrt((Ixx-Iyy)+4*Ixy^2);
            b=(lam-Ixx)/Ixy;                    %slope of division line
            d=Cy+b*(vx-Cx);                     %division line
            
            NRc=NRc+1;                          %number of cells increases by 1
            
            xt=find(vy>d);
            v=vx(xt)+(vy(xt)-1)*NVX;
            ctag(v)=NRc;
            csize(NRc)=length(xt);              %daughter cell has size of indices below division line
            csize(c)=csize(c)-length(xt);       %father cell shrinks by number of voxels to daughter cell
            
            %add intracellular parameters to this SUH 112220
            %8 intracellular concentrations- state variables
            % Update Statevars for cells by replicating cell c to end of
            % matrix
            new_statevars = nan(size(statevars,1)+8,1);
            %Put Extracellular state variables first: 4*NXY*NVY
            new_statevars(1:3*NVX*NVY,:) = statevars(1:3*NVX*NVY,:);  
            %Add variables that match the parent cell.
            intra=statevars(NVX*NVY*3+1:end,1);
            intra1=reshape(intra,NRc-1,8);
            intra1(end+1,:)=intra1(c,:);
            intra2=reshape(intra1,8*NRc,1);
            
            new_statevars(NVX*NVY*3+1:end, :) = intra2;
            statevars = new_statevars;
  
            % Update parammatrix for cells by replicating parameters for cell c to end of
            % parammatrix
            oldparammat=parammatrix;
            newparammat=oldparammat;
            newparammat(:,end+1)=oldparammat(:,c);
            parammatrix=newparammat;
            
            % Pass over proliferation rate to new cell
            phenoparams.pdivideE(end+1) = phenoparams.pdivideE(c); % SUH 3/24/20- OG was 0.003
            phenoparams.pdivideM(end+1) = phenoparams.pdivideM(c); %phenoparams.pdivideE(end)/3;
            cellmod.pop_index(end+1)=cellmod.pop_index(c);
            %pdivs(end+1)=pdivs(c); % Pass over proliferation rate to new cell
        end
    end
end


end
