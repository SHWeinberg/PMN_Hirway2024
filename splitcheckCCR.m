function split = splitcheckCCR(ctag,  csize,  xt, ttag, def)
NVX = def.NVX;
NV = def.NV;
greys = zeros(csize(ttag),1);
nbs(1)=xt-1+NVX; nbs(2)=xt+NVX; nbs(3)=xt+1+NVX;
nbs(8)=xt-1;                    nbs(4)=xt+1;
nbs(7)=xt-1-NVX; nbs(6)=xt-NVX; nbs(5)=xt+1-NVX;

prev = ctag(nbs(8)); in = 0;
for n = 1:8 % check Moore neighborhood to find target cell locations
    curr = ctag(nbs(n));
    if ((prev ~= ttag) && (curr == ttag)) % find the edges of the cell
        in = in + 1;
    end
    prev = curr;
end

split = false;
if in > 1 % the cell is at least one pixel in size

    % CONNECTED COMPONENT ALGORITHM Rene-style (CCR)
    % connected checking "label":
    % 0:  blue;   neighbors of retracted element
    % 1:  white;  undiscovered
    % 2:  grey;   discovered but not finished processing
    % 3:  black;  finished processing

    for v = 1:NV
        CCAlabels(v) = 1; % every pixel is an undiscovered location
    end

    CCAlabels(xt) = 3; % the current cell is the only processed location

    nrblue = -1;
    % cycle through Moore neighbors and find the current cell positions
    for n = 1:8
        nb = nbs(n);
        if ctag(nb) == ttag % the neighbor is self
            CCAlabels(nb) = 0; % the neighbor is a retracted element
            nrblue = nrblue + 1; % more than one neigbor must be self
            startnb = nb;
        end
    end
    CCAlabels(startnb) = 2; nrgrey = 1; greys(1) = startnb;
    % more than one neighbor is self and the Moore neighborhood needs processing
    while nrgrey && nrblue
        nrgrey0 = nrgrey;
        % make neighbors of discovered greys grey
        for i=1:nrgrey0
            g = greys(i);
            nbsg(1)=g-1+NVX; nbsg(2)=g+NVX; nbsg(3)=g+1+NVX;
            nbsg(8)=g-1;                    nbsg(4)=g+1;
            nbsg(7)=g-1-NVX; nbsg(6)=g-NVX; nbsg(5)=g+1-NVX;
            for n=1:8
                nb = nbsg(n);
                if (ctag(nb) == ttag) && (CCAlabels(nb) < 2)
                    if (CCAlabels(nb)==0)
                        nrblue = nrblue - 1;
                    end
                    CCAlabels(nb) = 2;
                    nrgrey = nrgrey + 1;
                    greys(nrgrey) = nb;
                end
            end
        end

        % make processed greys black
        for i = 1:nrgrey0
            g = greys(i);
            CCAlabels(g)=3;
            greys(i) = greys(nrgrey); nrgrey = nrgrey - 1;
        end
    end
    if nrblue
        split = true;
    end

end
end
