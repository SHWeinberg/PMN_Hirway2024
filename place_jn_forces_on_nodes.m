function [jfx, jfy]=place_jn_forces_on_nodes(ctag,celljn,NRc,jx,jy,def)
NV=def.NV;
NVX=def.NVX;

jfx=zeros(NV,1);
jfy=zeros(NV,1);
for c=1:NRc
    cellind=find(celljn==c);
    nind=length(cellind);
    for i=1:nind
        xs=cellind(i);
        stag=ctag(xs);
        nbs=[xs+NVX; xs+1; xs-NVX; xs-1];
        for n=1:4
            nb=nbs(n);
            nbtag=ctag(nb);
            if (nbtag~=stag) && (nbtag>0)
                p=nbtag+(stag-1)*NRc;
                jfx(xs)=jfx(xs)+jx(p);
                jfy(xs)=jfy(xs)+jy(p);
            end
        end
%         j(xs)=sqrt(jfx(xs)^2+jfy(xs)^2);
    end
end

end
