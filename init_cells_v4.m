function [NRc, ctag, csize, statevars, phenotype] = ...
    init_cells_v4(def,cellmod)
% With Initial matrix and ECM provided
NV = def.NV;
NVY = def.NVY;
NVX = def.NVX;
TARGETVOLUME = def.TARGETVOLUME;

icstate = cellmod.icstate;  % initial state of each cells
nstate = cellmod.nstate;  % number of state variables per cell
imatstate=cellmod.imatstate; % initial state of extracellullar concentrations
spread=cellmod.spread;
% ECAD_DENSITY = def.ECAD_DENSITY;

% ctag=zeros(NV,1);
% % ecadbzeros(NV,1);
% % ecadr=zeros(NV,1);
% NRc = 0; % plates cells
% for vy = 1:NVY
%     for vx = 1:NVX
%         v = vx + (vy - 1) * NVX;
%         if ((vx > 1) && (vx < NVX) && (vy > 1) && (vy < NVY)) % exclude outer rim
%             if (rand < (spread/TARGETVOLUME)) %smaller value makes it more sparse- SUH 112619 OG=0.25
%                 NRc = NRc + 1;
%                 ctag(v) = NRc;
%                 %                 ecadr(v) = ECAD_DENSITY;
%             end
%         end
%     end
% end
% 
% csize=zeros(NRc,1);
% for c = 1:NV
%     if (ctag(c) > 0)
%         csize(ctag(c)) = csize(ctag(c)) + 1;
%     end
% end
% 
% statevars = nan(3*NVX*NVY+NRc*nstate,1);
% matrvars=nan(NVX*NVY,3);
% matrvars(:,1)=imatstate(2); % solE
% matrvars(:,2)=imatstate(3)+imatstate(4); %totalE = asmE + ecmT
% matrvars(:,3)=imatstate(1)+imatstate(4); % totalT = solT + ecmT
% %matrvars(:,4)=imatstate(4);
% matrvars=reshape(matrvars,3*NVX*NVY,1);
% 
% cellvars=nan(NRc,nstate);
% for i = 1:NRc
%    cellvars(i,:) = icstate; 
% end
% cellvars=reshape(cellvars,NRc*nstate,1);
% statevars=[matrvars;cellvars];

% Providing intracellular and extracellular concentrations and initial cell seeding SUH 042921 
ctag=cellmod.initmatrix;
NRc=max(max(ctag));
csize=zeros(NRc,1);
for cell=1:NRc
    for c = 1:NV
        if (ctag(c) == cell)
            csize(cell) = csize(cell) + 1;
        end
    end
end

cellvars=cellmod.cellvars;
matrvars=cellmod.matrvars;
statevars=[matrvars;cellvars];

phenotype = calc_phenotype(statevars, cellmod.pheno,NRc);


end
