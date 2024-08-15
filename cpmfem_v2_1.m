function cpmfem_v2_1(filename, def, cellmod, phenoparams, cell_data)
warning('off','MATLAB:singularMatrix'), warning('off','MATLAB:rankDeficientMatrix');

%LOAD PARAMETER FILE
if nargin == 0
    def=init_vars();
end

%RNG SEED
rng(def.scurr);

% INITIALIZE
[NRc,ctag,csize, statevars, phenotype] = init_cells_v4(def, cellmod); % changed to v4 070721
[targetvols, pdivs, jccs, jcms] = calc_pheno_params_v3(phenotype, phenoparams);
[fx,fy,ux,uy] = set_forces(def);
[restrictx,restricty] = set_restrictions(def);
parammatrix=cellmod.parammatrix;

%LOCAL STIFFNESS MATRIX
klocal = set_klocal(def);

% GLOBAL STIFFNESS MATRIX
[kcol,kval] = assembly(klocal, def);
[dofpos, nrrdof] = arrange_dofpos(restrictx, restricty, def);
[kcol, kval] = reduce_K(kcol, kval, dofpos, nrrdof, def);

tic
for incr=1:def.NRINC
    disp(['Incr: ' num2str(incr)])
    %CPM: cell migration/proliferation
    [ctag, csize] = CPM_moves_v2(ctag, ux, uy, csize, def, targetvols, jccs, jcms,phenotype,phenoparams,cellmod);
    %[ctag,NRc,csize,statevars]=cell_proliferation_v2(ctag,NRc,csize,...
       % def,statevars,pdivs, targetvols); % Cell proliferation- properties
       % for new cells updated here- v2
  [ctag,NRc,csize,statevars,parammatrix,phenoparams,cellmod]=cell_proliferation_v3(ctag,NRc,csize,...
        def,statevars,pdivs, targetvols,parammatrix,cellmod,phenoparams); % Cell proliferation- properties for new cells updated here- v3 SUH 072221
    %FEA: calculate traction forces
    if def.fmaflag  % multicellular FMA
        [atag,cellmap,conn_mat_i,celljn] = connectedcomponent(ctag,NRc,def);
        [fx,fy]=cell_forces(atag,NRc,fx,fy,def);
        conn_mat = conn_mat_i;
    else  % single cell FMA
        [~,~,conn_mat_i,~] = connectedcomponent(ctag,NRc,def);
        atag=ctag;
        conn_mat=eye(NRc);
        cellmap=[1:NRc].*conn_mat;
        celljn=zeros(def.NV,1);
        [fx,fy]=cell_forces(ctag,NRc,fx,fy,def);
    end
    f = place_node_forces_in_f(fx,fy,restrictx,restricty,nrrdof,def);
    
    %FEA: calculate junction forces
    [jx,jy,jxnet,jynet,fxnet,fynet] = junction_forces_v2b(ctag, atag, NRc, ...
        conn_mat, cellmap, def);
    [jfx,jfy] = place_jn_forces_on_nodes(ctag,celljn,NRc,jx,jy,def);
    
    %FEA: calculate substrate strain
    u = set_disp_of_prev_incr(ux,uy,restrictx,restricty,nrrdof,def);
    u = solvePCG(kcol,kval,u,f,nrrdof,def);
    [ux,uy] = disp_to_nodes(u,restrictx, restricty, def);
    pstrain=get_pstrain(ux,uy,def);
    
    cellmod.params.mask=ctag; % Ctag mask passed into params- SUH112320
    %cellparameters=ones(NRc,1);
    %all_params=tian_model_baseline_params_v2(cellmod.tmaxval,cellmod.Jhalf,cellmod.extrascale1,cellmod.extrascale2,cellmod.intrazeb1,cellmod.intrasnail2,cellmod.pZEB,cellmod.pR200,cellmod.kdae_scale,NRc);
    %cellmod.params.model_params= all_params; - older v2
    % Can change individual cell parameters here- SUH 112320
        % New version of model baseline params-v3 with parameter heterogeneity
    % Add parameter for new cells- v3
    %Changed to v4 to fit parameter order
   all_params=tian_model_baseline_params_v4(parammatrix,cellmod.tmaxval,cellmod.Jhalf,cellmod.extrascale1,cellmod.extrascale2,cellmod.intrazeb1,cellmod.intrasnail2,cellmod.pZEB,cellmod.pR200,cellmod.kdae_scale,cellmod.ke_scale,cellmod.pop_index,NRc);
    cellmod.params.model_params= all_params; 
    
    % Cell ODE: update statevars/phenotype
    cellmod.params.cell_params = calc_cell_params(cell_data, conn_mat, jx, jy);  
    statevars = update_cellmarkers_v2(incr, statevars, cellmod.params,NRc);
    phenotype = calc_phenotype(statevars, cellmod.pheno,NRc);
    
   % save('phenotype1023.mat','phenotype');
   % save('phenoparams1023.mat','phenoparams');
    [targetvols, pdivs, jccs, jcms] = calc_pheno_params_v3(phenotype, phenoparams);

    % store variables
    sim_data.ctag(:,incr)=ctag;
    sim_data.csize{1,incr}=csize;
    sim_data.NRc(1,incr)=NRc;
    sim_data.conn_mat{1,incr}=sparse(conn_mat_i);
    sim_data.pstrain(:,incr)=pstrain(:,1);
    sim_data.pstrainx(:,incr)=pstrain(:,2);
    sim_data.pstrainy(:,incr)=pstrain(:,3);
    sim_data.jx{1,incr}=sparse(jx);
    sim_data.jy{1,incr}=sparse(jy);
    sim_data.jfx(:,incr)=jfx;
    sim_data.jfy(:,incr)=jfy;
    sim_data.jxnet{1,incr}=jxnet;
    sim_data.jynet{1,incr}=jynet;
    sim_data.fx(:,incr)=fx;
    sim_data.fy(:,incr)=fy;
    sim_data.fxnet{1,incr}=fxnet;
    sim_data.fynet{1,incr}=fynet;
    sim_data.statevars{1,incr} = statevars;
    sim_data.phenotype{1,incr} = phenotype;
    sim_data.parammatrix{1,incr}=parammatrix; % Stored parameter matrix to keep track of cell types
    sim_data.pop_index{1,incr}=cellmod.pop_index; % Stored cell types 
    
    
end
sim_data.Jcc=cellmod.JccTab; % Stored Jcc
sim_data.Jca=cellmod.JcaTab; % Stored Jca
sim_data.Jcm=cellmod.JcmTab; % Stored Jcm
    
sim_data.jfx = sparse(sim_data.jfx);
sim_data.jfy = sparse(sim_data.jfy);

write_data_v2(filename, sim_data,def, cellmod, phenoparams, cell_data);

end

