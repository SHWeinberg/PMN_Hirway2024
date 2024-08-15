function pstrain=get_pstrain(ux,uy,def)
NV=def.NV;
pstrain=zeros(NV,4);
estrains=[];
v1=[]; v2=[];
for v=1:NV
    estrains = get_estrains(ux, uy, v, estrains, def);
    L1=.0; L2=.0;
    [L1,L2,v1,v2] = get_princs(estrains,L1,L2,v1,v2,1);
    if L1>L2
        pstrain(v,1)=1E6*L1;
        pstrain(v,2)=(1E3*v1(1));
        pstrain(v,3)=(1E3*v1(2));
        pstrain(v,4)=1E6*L2;
    else
        pstrain(v,1)=1E6*L2;
        pstrain(v,2)=(1E3*v2(1));
        pstrain(v,3)=(1E3*v2(2));
        pstrain(v,4)=1E6*L1;
    end
end

end
