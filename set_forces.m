function [fx, fy, ux, uy] = set_forces(def)
NN = def.NN;
NNY = def.NNY;
NNX = def.NNX;
FORCE = def.FORCE;

a = (0.0/6.0) * 3.1416;

ux(NN, 1) = 0;
uy(NN, 1) = 0;

fx(NN, 1) = 0;
fy = fx;

for ny = 1:NNY
    for nx = 1:NNX
        n = nx + (ny-1)*NNX;

        % lower plate (iy==0) loading
        if ny == 1
            fx(n) = fx(n) + sin(a) * cos(a) * FORCE;
            fy(n) = fy(n) - cos(a) * cos(a) * FORCE;
        end
        % upper plate (iy==NNY-1) loading
        if ny == NNY
            fx(n) = fx(n) - sin(a) * cos(a) * FORCE;
            fy(n) = fy(n) + cos(a) * cos(a) * FORCE;
        end
        % left plate (ix==0) loading
        if nx == 1
            fx(n) = fx(n) - sin(a) * sin(a) * FORCE;
            fy(n) = fy(n) + sin(a) * cos(a) * FORCE;
        end
        % right plate (ix==NNX-1) loading
        if nx == NNX

            fx(n) = fx(n) + sin(a)*sin(a)*FORCE;
            fy(n) = fy(n) - sin(a)*cos(a)*FORCE;
        end
    end
end

for ny = 1:NNY
    for nx = 1:NNX
        n = nx + (ny-1)*NNX;
        % for loading on the side of a plate, forces are lower
        if (((nx == 1) || (nx == NNX)) && ((ny == 1) || (ny == NNY)))

            fx(n) = .5 * fx(n);
            fy(n) = .5 * fy(n);
        end
    end
end

end
