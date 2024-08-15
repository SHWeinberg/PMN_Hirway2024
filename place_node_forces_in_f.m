function f = place_node_forces_in_f(fx, fy, restrictx, restricty, nrrdof, def)
NN = def.NN;
f = zeros(nrrdof,1);
cnt = 0;
for n=1:NN
    if ~restrictx(n)
        cnt = cnt + 1;
        f(cnt) = fx(n);
    end
    if ~restricty(n)
        cnt = cnt + 1;
        f(cnt) = fy(n);
    end
end
end
