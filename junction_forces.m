function [jx,jy,jxnet,jynet,fxnet,fynet] = junction_forces(ctag, atag, NRc,...
                                            conn_mat, cellmap, def)
NNX=def.NNX;
NNY=def.NNY;
NVX=def.NVX;
VOXSIZE = def.VOXSIZE;
CELLFORCE=def.CELLFORCE;
NRj=NRc*(NRc-1)/2;  %#jnc=area of matrix - diagonal elements
fxnet=zeros(NRc,1); %net x-force for each cell
fynet=zeros(NRc,1); %net y-force for each cell
stop
% loop through cell clusters finding nodes
for a=1:NRc
    %get multicellular cluster nodes
    NRjxnn=0;
    for ny = 2:NNY - 1 % nodes 2 to NVY
        for nx = 2:NNX - 1 % nodes 2 to NVX
            n = nx + (ny - 1) * NNX; % nodes 302:90298
            anttag = 0;
            % look at the four nodes around the current node
            for vy = (ny - 1):ny % 1:300
                for vx = (nx - 1):nx % 1:300
                    v = vx + (vy - 1) * NVX;
                    if (atag(v) == a)
                        anttag = anttag + 1;
                    end
                end
            end
            if (anttag > 0)
                NRjxnn = NRjxnn + 1;
                jxnnodes(NRjxnn) = n;
            end
        end
    end
    
    %calculate force for each cell in cluster relative to cluster geometry
    cid=cellmap(a,cellmap(a,:)~=0);
    lim=length(cid);
    for i=1:lim
        %find each cell's nodes
        NRcelln=0;
        cellnodes=[];
        c=cid(i);
        % use inner nodes corresponding with voxels
        for ny = 2:NNY-1 % nodes 2 to NVY
            for nx = 2:NNX-1 % nodes 2 to NVX
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
        % calculate jxn force for each cell in cluster
        for p=1:NRcelln
            n=cellnodes(p);
            nx=rem(n-1,NNX)+1; ny=(n-nx)/NNY+1;
            for q=1:NRjxnn
                n2=jxnnodes(q);
                n2x=rem(n2-1,NNX)+1; n2y=(n2-n2x)/NNY+1;
                dny=(n2y-ny); % y distance between n and n2
                dnx=(n2x-nx); % x distance between n and n2
                fxnet(c)=fxnet(c)+dnx;
                fynet(c)=fynet(c)+dny;
            end
        end
    end
end
fxnet=fxnet*CELLFORCE*VOXSIZE;
fynet=fynet*CELLFORCE*VOXSIZE;

%JUNCTION MATRIX : RANK DEFICIENT
[c,m]=find(conn_mat);
n=c+(m-1)*NRc;
J=sparse(m,n,1,NRc,NRc^2);

%ASSUMPTION MATRIX OF REACTION FORCES
connu=triu(conn_mat,1); %upper triangle of connectivity matrix
[i,j]=find(connu);      %nonzero indices of junctions
p=j+(i-1)*NRc;          %value of 1
q=i+(j-1)*NRc;          %value of -1
k=1:length(p);
A = sparse(k,p,1,NRj, NRc^2) + sparse(k,q,1,NRj, NRc^2);

J=[J;A];
%T = [-fxnet -fynet; zeros(NRj,2)];
% x = J\T;
% jx = x(:,1);
% jy = x(:,2);

T = [-fxnet; zeros(NRj,1)];
jx = lsqminnorm(J,T);

T = [-fynet; zeros(NRj,1)];
jy = lsqminnorm(J,T);

jxnet=zeros(NRc,1);
jynet=zeros(NRc,1);
for i=1:NRc
    for j=1:NRc
        c=j+(i-1)*NRc;
        jxnet(i)=jxnet(i)+jx(c);
        jynet(i)=jynet(i)+jy(c);
    end
end

end
