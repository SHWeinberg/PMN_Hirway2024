function [restrictx, restricty] = set_restrictions(def)
NNX = def.NNX;
NNY = def.NNY;
NN = def.NN;

restrictx(NN, 1) = 0;
restricty(NN, 1) = 0;
for ny = 1:NNY
    for nx = 1:NNX
        if ((nx == 1) || (nx == NNX) || (ny == 1) || (ny == NNY))
            n = nx + (ny - 1) * NNX;
            restrictx(n) = true;
            restricty(n) = true;
        end
    end
end
end
