function [volumes, pdivs, jccs, jcms] = calc_pheno_params_test(phenotype, phenoparams)

% assumes a linear gradient between E and M cell phenotypes

% TARGETVOLUME function
mvol = (phenoparams.targetvolM - phenoparams.targetvolE);
volumes = phenoparams.targetvolE + mvol*phenotype;

% pdivide function
mdiv = (phenoparams.pdivideM - phenoparams.pdivideE);
pdivs = phenoparams.pdivideE + mdiv*phenotype;

% Jcc function
mjcc = (phenoparams.jccM - phenoparams.jccE);
jccs = phenoparams.jccE + mjcc*phenotype;

% Jcm function
mjcm = (phenoparams.jcmM - phenoparams.jcmE);
jcms = phenoparams.jcmE + mjcm*phenotype;

end