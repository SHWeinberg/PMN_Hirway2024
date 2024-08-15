function [dHcontact Hcontact Hcontactn] = calcdHcontact_v3(ctag, xt, ttag, stag, def,jccs, jcms,phenotype,phenoparams,cellmod)
NVX = def.NVX;

Hcontact=0; Hcontactn=0;
nbs(1)=xt-1+NVX; nbs(2)=xt+NVX; nbs(3)=xt+1+NVX;
nbs(8)=xt-1;                    nbs(4)=xt+1;
nbs(7)=xt-1-NVX; nbs(6)=xt-NVX; nbs(5)=xt+1-NVX;
%disp(['this is nbs: ' num2str(nbs)]);
if1=0;
else1=0;
for n=1:8
  %  disp(['This is N: ' num2str(n)])
    nbr = nbs(n);
   % disp(['this is nbr: ' num2str(nbr)]);
    nbtag = ctag(nbr);
   % disp(['this is nbtag: ' num2str(nbtag)]);
   
   
     if (stag==0 && nbtag~=0 && nbtag ~=ttag && ttag~=0)
       Hcontactn = Hcontactn + contactenergy_v3(stag, nbtag, jccs, jcms,phenotype,phenoparams,cellmod)+JcApheno(ttag,nbtag,phenotype,phenoparams,cellmod);
       % disp(['1.Contactenergy: ' num2str(contactenergy_v3(stag, nbtag, jccs, jcms,phenotype,phenoparams)+JcApheno(ttag,nbtag,phenotype,phenoparams))]);
       if1=if1+1;
      %  disp(['sTag & nbTag: ' num2str(stag) ', ' num2str(nbtag)])
      %  disp('___________')
    else
        Hcontactn = Hcontactn + contactenergy_v3(stag, nbtag, jccs, jcms,phenotype,phenoparams,cellmod);
      %  disp(['2.ContactEnergy: ' num2str(contactenergy_v3(stag, nbtag, jccs, jcms,phenotype,phenoparams))]);
      %  disp(['sTag & nbTag: ' num2str(stag) ', ' num2str(nbtag)])
        else1=else1+1;
      %     disp('___________')
    end
    Hcontact = Hcontact + contactenergy_v3(ttag, nbtag, jccs, jcms,phenotype,phenoparams,cellmod);
   %  disp(['3.ContactH: ' num2str(contactenergy_v3(ttag, nbtag, jccs, jcms,phenotype,phenoparams))]);
   % disp(['tTag & nbTag: ' num2str(ttag) ', ' num2str(nbtag)])
     %   disp('___________')
end
%disp(['if counter: ' num2str(if1)])
%disp(['else counter: ' num2str(else1)])
dHcontact = Hcontactn - Hcontact;

end
