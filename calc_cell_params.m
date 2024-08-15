function cell_params = calc_cell_params(cell_data, conn_mat, jx, jy)

method = cell_data.method;
TGFb_max = cell_data.TGFb_max;
TGFBval=cell_data.TGFBval;
[ncell,~] = size(conn_mat);
%New

%
switch method
    case 'TGFb_Conn'
        neighbor_max = cell_data.neighbor_max;
        num_neighbors = sum(sign(conn_mat));
        
        %New- SUH 112519
        jmag = sqrt(jx.^2 + jy.^2);
        jmag_mat = reshape(jmag, ncell, ncell);
        jmag_cell = sum(jmag_mat);   % magnitude of junctional forces on each cell
        %%
        for i = 1:ncell
           %disp(['Ncell: ' num2str(ncell)]);
            cell_params{i}.TGF0 = TGFBval; %max(0,TGFb_max - TGFb_max/neighbor_max*num_neighbors(i));
            % make 1 value, initiated in run cpm model file
            %SUH 112919
            %New- SUH 112519
            
            cell_params{i}.J=jmag_cell(i); %NEED TO ADD BACK IN SUH-120519
            %add cell_params{i}.Tmax;
            %             if num_neighbors(i)
            %                cell_params{i}.TGF0 = 0;
            %             else
            %                cell_params{i}.TGF0 = 3;
            %             end
        end
    case 'TGFb_Junc'
        jmag_max = cell_data.jmag_max;
        jmag = sqrt(jx.^2 + jy.^2);
        jmag_mat = reshape(jmag, ncell, ncell);
        jmag_cell = sum(jmag_mat);   % magnitude of junctional forces on each cell
        
        for i = 1:ncell
            cell_params{i}.TGF0 = max(0, TGFb_max - TGFb_max*jmag_cell(i)/jmag_max);
        end
end

end