function phenotype = calc_phenotype(statevars, pheno,NRc)
% [ncell,~] = size(statevars);

method = pheno.method;

switch method
    case 'NcadConc'
        phenotype = statevars(end-NRc+1:end)/pheno.ncad_max; 
       % disp(['ncad: ' num2str(pheno.ncad_max)]);
        %phenotype=phenotype*0; % change to 0 to make cell Epithelial-
        %Commented out to have phenotype dependence again SUH-112919
end


end