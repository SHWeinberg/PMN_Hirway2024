function [fx, fy] = cell_forces(ctag, NRc, fx, fy, def)
NVX=def.NVX;
NNX=def.NNX;
NNY=def.NNY;
CELLFORCE = def.CELLFORCE;
VOXSIZE = def.VOXSIZE;

%reset inner matrix to 0
for ny=2:NNY-1
    for nx=2:NNX-1
        n = nx + (ny-1) * NNX;
        fx(n) = 0;
        fy(n) = 0;
    end
end

for c=1:NRc
    NRcelln=0;
    % use inner nodes corresponding with voxels
    for ny = 2:NNY - 1 % nodes 2 to NVY
        for nx = 2:NNX - 1 % nodes 2 to NVX
            n = nx + (ny - 1) * NNX; % nodes 302:90298
            cnttag = 0;
            % look at the box around the current node
            for vy = (ny - 1):ny % 1:300
                for vx = (nx - 1):nx % 1:300
                    v = vx + (vy - 1) * NVX;
                    if (ctag(v) == c)
                        cnttag = cnttag + 1; % one of the box is self
                    end
                end
            end
            if (cnttag > 0)
                NRcelln = NRcelln + 1;
                cellnodes(NRcelln) = n;
            end
        end
    end
    for i=1:NRcelln
        n=cellnodes(i);
        nx=rem(n-1,NNX)+1; ny=(n-nx)/NNY+1;
        for j=1:NRcelln
            n2=cellnodes(j);
            n2x=rem(n2-1,NNX)+1; n2y=(n2-n2x)/NNY+1;
            dnx=n2x-nx; % x distance between n and n2
            dny=n2y-ny; % y distance between n and n2
            fx(n)=fx(n)+dnx;
            fy(n)=fy(n)+dny;
        end
    end
end
fx=fx*CELLFORCE*VOXSIZE;
fy=fy*CELLFORCE*VOXSIZE;
end
