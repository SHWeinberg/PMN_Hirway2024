function [L1,L2,v1,v2] = get_princs( str,  L1,  L2,  v1,  v2,  strain)

xx=str(1); yy=str(2);
if strain
    xy = 0.5*str(3);
else
    xy = str(3);
end
% mat = xx xy
%       xy yy


if xy == 0 % no shear strain
    L1 = xx; v1(1)=1; v1(2)=0;
    L2 = yy; v2(1)=0; v2(2)=1;
else

    T = xx + yy; % trace
    D = xx*yy - xy*xy; % determinant

    T2D = T*T/4-D;

    if T2D <= 0 % can occur if strain very close to isotropic

        L1 = xx; v1(1)=1; v1(2)=0;
        L2 = yy; v2(1)=0; v2(2)=1;

    else

        sqT2D = sqrt(T2D);
        L1 = T/2+sqT2D;
        L2 = T/2-sqT2D;

        % eigenvector v must satisfy: mat*v = L*v
        % xx*v[0]+xy*v[1] = L*v[0] -> if v[0]=1, v[1]=(L-xx)/xy (=Q)
        % ||v|| = sqrt(1+Q*Q) (=R)   -> v = [1/R;Q/R]
        Q=(L1-xx)/xy; R=sqrt(1+Q*Q); v1(1)=1/R; v1(2)=Q/R;
        Q=(L2-xx)/xy; R=sqrt(1+Q*Q); v2(1)=1/R; v2(2)=Q/R;
    end
end
end
