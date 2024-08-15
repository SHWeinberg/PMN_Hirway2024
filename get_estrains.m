function estrains =  get_estrains(ux, uy, e, estrains,def)
NVX = def.NVX;
NVY = def.NVY;
NNX = def.NNX;
B = zeros(3,8);
B = set_matrix_B(B, 0, 0,def);
[vx,vy]=ind2sub([NVX NVY],e);
% vx=rem(e,NVX) + 1; vy=(e-vx)/NVX + 1;
% determine corner node numbers of this voxel
n00 = (vx  ) + (vy - 1)*NNX;
n10 = (vx+1) + (vy - 1)*NNX;
n11 = (vx+1) + (vy)*NNX;
n01 = (vx  ) + (vy)*NNX;

u(1) = ux(n00);
u(2) = uy(n00);
u(3) = ux(n10);
u(4) = uy(n10);
u(5) = ux(n11);
u(6) = uy(n11);
u(7) = ux(n01);
u(8) = uy(n01);

for i=1:3
    estrains(i) = 0;
end
% strain displacement relation
for i=1:3
    for j=1:8
        estrains(i) = estrains(i) + B(i,j) * u(j);
    end
end
end
