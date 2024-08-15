function klocal = set_klocal(def)
VOXSIZE = def.VOXSIZE;

% two-point Gaussian integration
% local coordinates of the integration points (1/sqrt(3) = 0.57735027)
intgrx = [-.57735,  .57735,  .57735, -.57735];
intgry = [-.57735, -.57735,  .57735,  .57735];

D = zeros(3,3);		% material matrix
B = zeros(3,8);		% strain displacement matrix
Bt = zeros(8,3);	% transpose of strain displacement matrix
BD = zeros(8,3);	% Bt * D
BDB = zeros(8,8);	% BD * B

% node positions in local coordinate system
% double nx[4] = {-1,  1,  1, -1};
% double ny[4] = {-1, -1,  1,  1};

% set matrix k to zeros
for m=1:8
    for n=1:8
        klocal(m,n) = 0.0;
    end
end

% calculate stifness matrix of the material (linear elastic isotropic)
D = material_matrix(D,def);

% Determine local stiffness matrix
% by integration. Implemented as summation over all integration points
for i=1:4      % for all integration points
    % calculate matrix B in intgr pnt i
    B = set_matrix_B(B, intgrx(i), intgry(i),def);
    
    % Bt is the transpose of B
    for m=1:8
        for n=1:3
            Bt(m,n) = B(n,m);
        end
    end
    % BD  =  Bt * D
    for m=1:8
        for n=1:3
            BD(m,n) = 0;
            for j=1:3
                BD(m,n) = BD(m,n) + Bt(m,j) * D(j,n);
            end
        end
    end
    % BDB  =  BD * B
    for m=1:8
        for n=1:8
            BDB(m,n) = 0;
            for j=1:3
                BDB(m,n) = BDB(m,n) + BD(m,j) * B(j,n);
            end
        end
    end
    
    % integration over the volume. This leads to adding to local
    % stifness matrix for each integration point
    % dV = dx*dy*dz = det(J) * dr*ds*dt
    % for cubic voxel elements the volume represented by one
    % integration point is equal to dV = (.5*VOXSIZE)^3
    dV = .25 * VOXSIZE * VOXSIZE;
    for m=1:8
        for n=1:8
            klocal(m,n) = klocal(m,n) + BDB(m,n) * dV;
        end
    end
end % endfor all integration points
end
