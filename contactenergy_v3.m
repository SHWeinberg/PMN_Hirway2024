function J = contactenergy_v3(tag1, tag2, jccs, jcms,phenotype,phenoparams,cellmod)
%   JCM = def.JCM;
%   JCC = def.JCC;

  J = 0; tag = max(tag1, tag2);
 % disp(['this is tag:' num2str(tag)]) %SUH 10/11/19
  % disp(['this is tag1:' num2str(tag1)]) %SUH 10/11/19
   % disp(['this is tag2:' num2str(tag2)]) %SUH 10/11/19
  if tag1 ~= tag2
    if (tag1 == 0) || (tag2 == 0)
     % J = jcms(tag);
     J=Jcmpheno(tag,phenotype,phenoparams,cellmod);
    % disp('Jcm')
    else
      %jccs case
      %J = jccs(tag); %SUH-10/23
      %disp(['Tag1 & Tag2: ' num2str(tag1) ', ' num2str(tag2)])
      %disp(['Phenotype: ' num2str(phenotype')])
      J=Jccpheno(tag1,tag2,phenotype,phenoparams,cellmod);
    %  disp('Jcc')
    end
  end
