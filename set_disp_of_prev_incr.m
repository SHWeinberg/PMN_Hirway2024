function u = set_disp_of_prev_incr(ux, uy, restrictx, restricty, nrrdof, def)
NN = def.NN;

u=zeros(nrrdof,1);
cnt = 0;
for n=1:NN
    if ~restrictx(n)
        cnt = cnt + 1;
        u(cnt) = ux(n);
    end
    if ~restricty(n)
        cnt = cnt + 1;
        u(cnt) = uy(n);
    end
end
end
