function [ux, uy] = disp_to_nodes(u, restrictx, restricty, def)
NN = def.NN;

ux=zeros(NN,1);
uy=zeros(NN,1);
cnt = 0;
for n = 1:NN
    if ~restrictx(n)
        cnt = cnt + 1;
        ux(n) = u(cnt);
    end
    if ~restricty(n)
        cnt = cnt + 1;
        uy(n) = u(cnt);
    end
end
end
