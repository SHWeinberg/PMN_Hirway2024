function [kcol,kval]=reduce_K(kcol,kval,dofpos,nrrdof,def)
NDOF = def.NDOF;

for ro = 1:NDOF
    rn = dofpos(ro);
    if rn > -1 % if this row is not to be removed
        lim = kcol(ro, 1);
        %         lim = 10*(ro-1)+kcol(10*(ro-1)+1);
        for a = 2:lim
            %         for a=(10*(ro-1)+2):lim
            co = kcol(ro, a); % old column index
            cn = dofpos(co);  % new column index
            kcol(ro, a) = cn; % give new column index (some get -1)
        end
        % remove columns with -1 index
        shift = 0;
        for a = 2:lim
            %         for a=(10*(ro-1)+2):lim
            kcol(ro, a-shift) = kcol(ro, a);
            kval(ro, a-shift) = kval(ro, a);
            if kcol(ro, a) == -1
                shift = shift + 1;
            end
        end
        %         kcol(10*(ro-1)+1) = kcol(10*(ro-1)+1)-shift;
        kcol(ro, lim-shift+1:lim)=0;
        kcol(ro, 1) = kcol(ro, 1) - shift;
        % shift row itself
        for a=1:10
            kcol(rn, a) = kcol(ro, a);
            kval(rn, a) = kval(ro, a);
            %             kcol(10*(rn-1)+a) = kcol(10*(ro-1)+a);
            %             kval(10*(rn-1)+a) = kval(10*(ro-1)+a);
        end
    end
end
%clear shifted columns
kcol(nrrdof+1:end,:) = [];
kval(nrrdof+1:end,:) = [];

end
