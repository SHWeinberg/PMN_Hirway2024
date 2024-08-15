function dHstrain = calcdHstrain(ux, uy, xt, xs, pick, ttag, stag,def)
SQ05 = def.SQ05;
YOUNGS = def.YOUNGS;
STIFFENINGSTIFF = def.STIFFENINGSTIFF;
v1 = []; v2 = [];
estrains = [];
dHstrain = 0;

% unitvectors for move: vm
q = SQ05;
vmx(1)=-q; vmx(2)= 0; vmx(3)= q;
vmx(8)=-1;            vmx(4)= 1;
vmx(7)=-q; vmx(6)= 0; vmx(5)= q;

vmy(1)= q; vmy(2)= 1; vmy(3)= q;
vmy(8)= 0;            vmy(4)= 0;
vmy(7)=-q; vmy(6)=-1; vmy(5)=-q;

% eigenvector = strain orientation
vm(1) = vmx(pick);
vm(2) = vmy(pick);


if stag % expansion

    estrains = get_estrains(ux, uy, xt, estrains, def);
    L1=.0; L2=.0; [L1,L2,v1,v2] = get_princs(estrains,L1,L2,v1,v2,1);

    % inproducts of move vector with pr. strain vectors
    vmv1 = vm(1)*v1(1) + vm(2)*v1(2);
    vmv2 = vm(1)*v2(1) + vm(2)*v2(2);

    %dHstrain -= sige(L1)*vmv1*vmv1 + sige(L2)*vmv2*vmv2;
    %dHstrain -= sige(YOUNGS);

    E1 = YOUNGS;
    if L1 > 0
        E1 = E1*(1+L1/STIFFENINGSTIFF);
    end
    E2 = YOUNGS;
    if L2 > 0
        E2 = E2 * (1+L2/STIFFENINGSTIFF);
    end
    dHstrain = dHstrain - sige(E1,def)*vmv1*vmv1 + sige(E2,def)*vmv2*vmv2;
end

if ttag % retraction

    estrains = get_estrains(ux, uy, xs, estrains, def);
    L1=.0; L2=.0; [L1,L2,v1,v2] = get_princs(estrains,L1,L2,v1,v2,1);

    % inproducts of move vector with pr. strain vectors
    vmv1 = vm(1)*v1(1) + vm(2)*v1(2);
    vmv2 = vm(1)*v2(1) + vm(2)*v2(2);

    %dHstrain += sige(L1)*vmv1*vmv1 + sige(L2)*vmv2*vmv2;
    %dHstrain += sige(YOUNGS);

    E1 = YOUNGS;
    if L1>0
        E1 = E1*(1+L1/STIFFENINGSTIFF);
    end

    E2 = YOUNGS;
    if L2>0
        E2 = E2*(1+L2/STIFFENINGSTIFF);
    end
    dHstrain = dHstrain + sige(E1,def)*vmv1*vmv1 + sige(E2,def)*vmv2*vmv2;
end

end
