function new_statevars = update_cellmarkers_v2(mcs, old_statevars, params,NRc)
% to match emt_tian_2013_v3 file

[NRc] = NRc;%size(old_statevars);
disp(['Cells: ' num2str(NRc)]);
%disp(['Statvars: ' num2str(size(old_statevars))]);


method = params.method;
dt = params.dt;  % time step, min
odefun = params.odefun;  % cell ODE function
ti = (mcs-1)*params.mcs_to_time;  % time, min
%model_params = params.model_params;  % structure w/model parameters (non-cell dependent)
par=struct;
par.parametermatrix=params.model_params.parametermatrix;
par.extra=params.model_params.extra;
par.mask=params.mask;
par.nx=size(params.mask,1);
par.cells=NRc;


par.scale=2;
par.conv=1;


cellTgfb=[];
cellJ=[];
for c=1:NRc % Create vectors for cell parameters
    cellTgfb(c)=params.cell_params{c}.TGF0; %Exotgfb
    cellJ(c)=params.cell_params{c}.J;
end
par.cellTgfb=cellTgfb';
par.cellJ=cellJ';
%par.exo=params.cell_
new_statevars = nan(size(old_statevars));
%for i = 1:NRc %get rid of, make it go through entire grid and cells
   % cell_params = params.cell_params{i};
    switch method
        case 'SimpleEuler'  % if dt = mcs_to_time
            new_statevars = old_statevars' + dt*odefun(ti, old_statevars, par);
        case 'Euler'
            neuler = params.mcs_to_time/dt;  % should be divisible
            statevar = old_statevars(i,:);
            statevar = statevar(:);
            for j = 1:neuler-1
                statevar = statevar + dt*odefun(ti+(j-1)*dt, statevar, model_params, cell_params);
            end
            new_statevars(i,:) = statevar;

        case 'MatODE'
            statevar = old_statevars;
            %statevar = statevar(:);
            fun = @(t,x) odefun(t, x, par);
            [~,X] = ode45(fun,[ti ti+params.mcs_to_time],statevar);
            X1=X(end,:)';
            new_statevars = X1;
    end