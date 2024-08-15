function cpmfem(varargin)
warning('off','MATLAB:singularMatrix'), warning('off','MATLAB:rankDeficientMatrix');

%LOAD PARAMETER FILE
if nargin>0
    for n=1:2:nargin
        switch varargin{n}
            case 'def'
                def=varargin{n+1};
        end
    end
else
    def=init_vars();
end

%RNG SEED
rng(def.scurr);

% INITIALIZE
[NRc,ctag,csize] = init_cells(def);
[fx,fy,ux,uy] = set_forces(def);
[restrictx,restricty] = set_restrictions(def);

%LOCAL STIFFNESS MATRIX
klocal = set_klocal(def);

% GLOBAL STIFFNESS MATRIX
[kcol,kval] = assembly(klocal, def);
[dofpos, nrrdof] = arrange_dofpos(restrictx, restricty, def);
[kcol, kval] = reduce_K(kcol, kval, dofpos, nrrdof, def);

tic
for incr=1:def.NRINC
    %FEA: calculate traction forces
    if def.fmaflag
        [atag,cellmap,conn_mat_i,celljn] = connectedcomponent(ctag,NRc,def);
        [fx,fy]=cell_forces(atag,NRc,fx,fy,def);
        conn_mat = conn_mat_i;
    else
        [~,~,conn_mat_i,~] = connectedcomponent(ctag,NRc,def);

        atag=ctag;
        conn_mat=eye(NRc);
        cellmap=[1:NRc].*conn_mat;
        celljn=zeros(def.NV,1);
        [fx,fy]=cell_forces(ctag,NRc,fx,fy,def);
    end
    f = place_node_forces_in_f(fx,fy,restrictx,restricty,nrrdof,def);
    
    %FEA: calculate junction forces
    [jx,jy,jxnet,jynet,fxnet,fynet] = junction_forces(ctag, atag, NRc, ...
        conn_mat, cellmap, def);
    [jfx,jfy] = place_jn_forces_on_nodes(ctag,celljn,NRc,jx,jy,def);
    
    %FEA: calculate substrate strain
    u = set_disp_of_prev_incr(ux,uy,restrictx,restricty,nrrdof,def);
    u = solvePCG(kcol,kval,u,f,nrrdof,def);
    [ux,uy] = disp_to_nodes(u,restrictx, restricty, def);
    pstrain=get_pstrain(ux,uy,def);
    
    %CPM: cell migration
    [ctag, csize] = CPM_moves(ctag, ux, uy, csize, def);
    [ctag,NRc,csize]=cell_proliferation(ctag,NRc,csize,def);
    
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
end

sim_data.jfx = sparse(sim_data.jfx);
sim_data.jfy = sparse(sim_data.jfy);

write_data(sim_data,def);

end

