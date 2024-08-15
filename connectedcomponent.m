function [atag,cellmap,conn_mat,celljn] = connectedcomponent(ctag, NRc, def)
NVX = def.NVX;
NV = def.NV;
cellind = find(ctag);
nind=length(cellind);
CCAlabels=ones(NV, 1);
atag=zeros(NV, 1);
NRa=NRc;

%flood-fill algorithm
for c=1:NRc
    nrblue = 0;
    for a = 1:nind
        xs = cellind(a);
        stag = ctag(xs);
        if (stag == c) && (CCAlabels(xs) < 2)
            nrblue = nrblue + 1;
            CCAlabels(xs) = 0;
            startnb = xs;
        end
    end
    nrgrey = 1; greys(1) = startnb; CCAlabels(startnb) = 2;
    while nrgrey && nrblue
        nrgrey0 = nrgrey;
        for i = 1:nrgrey0
            g = greys(i);
            nbsg(1)=g-1+NVX; nbsg(2)=g+NVX; nbsg(3)=g+1+NVX;
            nbsg(8)=g-1;                    nbsg(4)=g+1;
            nbsg(7)=g-1-NVX; nbsg(6)=g-NVX; nbsg(5)=g+1-NVX;
            for n = 1:8
                nb = nbsg(n);
                ttag = ctag(nb);
                if (ttag > 0) && (CCAlabels(nb) < 2)
                    if (CCAlabels(nb) == 0)
                        nrblue=nrblue-1;
                    end
                    CCAlabels(nb)=2;
                    nrgrey=nrgrey+1;
                    greys(nrgrey)=nb;
                end
            end
        end
        for i=1:nrgrey0
            g=greys(i);
            CCAlabels(g)=3;
            greys(i)=greys(nrgrey);
            nrgrey=nrgrey-1;
            atag(g)=c;
        end
    end
end

conn_mat=zeros(NRc,NRc);
celljn=zeros(NV,1);
cellmap=zeros(NRc,NRc);
for c=1:NRc
    for p=1:nind
        xs=cellind(p);
        stag=ctag(xs);
        if (stag==c)
            nbs(1)=xs-1+NVX; nbs(2)=xs+NVX; nbs(3)=xs+1+NVX;
            nbs(8)=xs-1;                    nbs(4)=xs+1;
            nbs(7)=xs-1-NVX; nbs(6)=xs-NVX; nbs(5)=xs+1-NVX;
            for n=1:8
                nb=nbs(n);
                nbtag=ctag(nb);
                if (nbtag>0) && (nbtag~=stag)
                    conn_mat(c,nbtag)=conn_mat(c,nbtag)+1;   %connectivity matrix
                end
            end
        end
        
        if (atag(xs)==c)
            cellmap(c,ctag(xs))=ctag(xs);   %cell mapping to clusters
        end
    end
end

CCAlabels=zeros(NV,1);
for c=1:NRc    
    lim=length(cellmap(c,cellmap(c,:)~=0));
    for p=1:lim
        stag=cellmap(c,p);
        for q=1:nind
            xs=cellind(q);
            ttag=ctag(xs);
            if (ttag==stag)
                nbs(1)=xs-1+NVX; nbs(2)=xs+NVX; nbs(3)=xs+1+NVX;
                nbs(8)=xs-1;                    nbs(4)=xs+1;
                nbs(7)=xs-1-NVX; nbs(6)=xs-NVX; nbs(5)=xs+1-NVX;
                for n=1:8
                    nb=nbs(n);
                    nbtag=ctag(nb);
                    if nbtag~=stag && nbtag>0
                        if CCAlabels(xs)<1
                            celljn(xs)=c;
                            CCAlabels(xs)=1;
                        end
                        if CCAlabels(nb)<1
                            celljn(nb)=c;
                            CCAlabels(nb)=1;
                        end
                    end
                end
            end
        end
    end
end


end
