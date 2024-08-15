function [kcol, kval] = assembly(klocal, def)
NDOF = def.NDOF;
NVY = def.NVY;
NVX = def.NVX;
NNX = def.NNX;

kcol(NDOF, 10) = 0; kcol(:, 1) = 1;
kval(NDOF, 10) = 0;
for vy=1:NVY
    for vx=1:NVX
        Ef = 1;%(vx+1)/(double)NVX+.5;
        
        % determine corner node numbers of this element
        n00 = (vx) + (vy - 1)*NNX;
        n10 = (vx+1) + (vy - 1)*NNX;
        n11 = (vx+1) + (vy)*NNX;
        n01 = (vx) + (vy)*NNX;
        
        topv(1) = 2*n00-1;
        topv(2) = 2*n00;
        topv(3) = 2*n10-1;
        topv(4) = 2*n10;
        topv(5) = 2*n11-1;
        topv(6) = 2*n11;
        topv(7) = 2*n01-1;
        topv(8) = 2*n01;
        
        % place klocal in K matrix
        for il=1:8 % go through rows in klocal
            for jl=1:8 % go through columns in klocal
                value = Ef * klocal(il,jl);
                ig = topv(il); % row in K
                jg = topv(jl); % column in K
                if jg == ig % if on diagonal
                    kval(ig, 1) = kval(ig, 1) + value;
                end
                if jg > ig % if right of diagonal
                    lim = kcol(ig, 1);
                    % check if there was already a nonzero on K(ig,jg)
                    alreadynonzero = false;
                    for a = 2:lim
                        if kcol(ig, a) == jg % if storage for jg-th column of K
                            alreadynonzero = true;
                            b = a;
                        end
                    end
                    if alreadynonzero % if already a nonzero on K(ig,jg)
                        kval(ig, b) = kval(ig, b) + value; % add klocal(il,jl)
                    else % if nothing on K(ig,jg)
                        b = lim + 1;
                        kcol(ig, b) = jg; % make storage for jg-th column of K
                        kval(ig, b) = value; % add klocal(il,jl)
                        kcol(ig, 1) = kcol(ig, 1) + 1;
                    end
                end %endfor go through klocal
            end
        end %endif relevant element
    end
end % endfor go though elements

end
