function [dofpos,cnt] = arrange_dofpos(restrictx, restricty, def)
NN = def.NN;
NDOF = def.NDOF;

dofpos(NDOF, 1) = 0;
cnt = 0;
for n=1:NN
    if (restrictx(n))
        dofpos(2*n-1) = -1;
    else
      cnt = cnt + 1;
      dofpos(2*n-1) = cnt;
    end
    if (restricty(n))
        dofpos(2*n) = -1;
    else
        cnt = cnt + 1;
        dofpos(2*n) = cnt;
    end
end
end
