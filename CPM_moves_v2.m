function [ctag, csize] = CPM_moves_v2(ctag, ux, uy,  csize, def, ...
    targetvols, jccs, jcms,phenotype,phenoparams,cellmod)
IMMOTILITY = def.IMMOTILITY;
NVX = def.NVX;
NVY = def.NVY;
NV = def.NV;
NRsteps = NV;
%   ECAD_DENSITY = def.ECAD_DENSITY;
% cellular potts model: one Monte Carlo step
for i = 1:NRsteps
    xt = ceil(rand * NV); % pick random element
    xtx = rem(xt - 1, NVX) + 1; xty = (xt - xtx)/NVY + 1;
    
    
    if ((xtx > 1) && (xtx < NVX) && (xty > 1) && (xty < NVY)) % exclude outer rim
        nbs(1)=xt-1+NVX; nbs(2)=xt+NVX; nbs(3)=xt+1+NVX;
        nbs(8)=xt-1;                    nbs(4)=xt+1;
        nbs(7)=xt-1-NVX; nbs(6)=xt-NVX; nbs(5)=xt+1-NVX;
        pick = ceil(rand * 8);
        xs = nbs(pick); % source pixel
        
        ttag = ctag(xt);
        stag = ctag(xs);
        
        go_on = 0;
        if (ttag ~= stag) % don't bother if no difference
            go_on = 1;
            if ttag % if a cell in xt (retracting)
                split = splitcheckCCR(ctag, csize, xt, ttag, def);
                if split
                    go_on = 0;
                end
                if csize(ttag) == 1 % cell cannot disappear (constraint may be removed)
                    go_on = 0;
                end
            end
        end
        
        if go_on
            dH = calcdH_v2(ctag, ux, uy, csize, xt, xs, pick, ttag, stag, def, ...
                targetvols, jccs, jcms,phenotype,phenoparams,cellmod);
            prob = exp(-IMMOTILITY * dH);
            if prob > rand()
                ctag(xt) = stag; % a move is made
                
                % cell in xt is retracting
                if ttag
                    csize(ttag) = csize(ttag) - 1;
                    
                end
                
                % cell in xs is extending into xt
                if stag
                    csize(stag) = csize(stag) + 1;
                    
                end
            end
        end
    end
end
