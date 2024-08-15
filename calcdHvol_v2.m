function dHvol = calcdHvol_v2(csize, ttag, stag, def, targetvols)
% TARGETVOLUME = def.TARGETVOLUME;
INELASTICITY = def.INELASTICITY;

dHvolA = 0; dHvolB = 0; 

if ttag % cell ttag retracts
    V = csize(ttag); 
    V0 = targetvols(ttag);
    eV = (V - V0)/V0; eVn = (V - 1 - V0)/V0;
    dHvolA = INELASTICITY * ((eVn * eVn) - (eV * eV));
end

if stag % cell stag expands
    V = csize(stag); 
    V0 = targetvols(stag);
    eV = (V - V0)/V0; eVn = (V + 1 - V0)/V0;
    dHvolB = INELASTICITY * ((eVn * eVn) - (eV * eV));
end

dHvol = dHvolA + dHvolB;

end
