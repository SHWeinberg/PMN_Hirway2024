function [xjunc_cent, yjunc_cent, jx_cent, jy_cent, ...
    xjunc_vox, yjunc_vox, jx_vox_norm, jy_vox_norm, jx_vox, jy_vox,...
    jx_mat, jy_mat, neighbors_mat, Aconn] ...
    = calc_junc_from_ctag(ctag_mat, nx, ny)

% warning('off','MATLAB:singularMatrix'), warning('off','MATLAB:rankDeficientMatrix');

icell = unique(ctag_mat(:));
icell = icell(icell>0);

% Ncell = length(icell);
Ncell = max(ctag_mat(:));

% calculate connectivity matrix
Aconn = zeros(Ncell);
for i = 1:length(icell)  % loop thru each cell
    [y,x] = find(ctag_mat==icell(i));  % voxels for cell i
    for p = 1:length(y)  % loop thru each voxel in cell i
        xi = x(p); yi = y(p);
        % loop through all 4 neighbors
        for k = 1:4
            switch k
                case 1
                    xn = xi-1; yn = yi;
                case 2
                    xn = xi+1; yn = yi;
                case 3
                    xn = xi; yn = yi-1;
                case 4
                    xn = xi; yn = yi+1;
            end
            possible_neighbor = ctag_mat(yn, xn);
            if possible_neighbor>0 && possible_neighbor~=icell(i)
                Aconn(icell(i),possible_neighbor) = 1;
                Aconn(possible_neighbor, icell(i)) = 1;
            end
            
        end
        
        
    end
end

% find clusters
G = graph(Aconn);
bins = conncomp(G);  % finds connected components from connectivity matrix
% bins is a 1 x Ncell vector, where each element is the cluster index for
% each cell
num_clusters = max(bins);  

atag_mat = zeros(size(ctag_mat));
fxnet = zeros(Ncell,1);  % net traction forces for each cell
fynet = zeros(Ncell,1);
cell_centroids = zeros(Ncell,2);

% vectors of traction forces for all nodes
traction_x = []; traction_y = [];  % positions
traction_u = []; traction_v = [];  % traction vectors

cluster_centroids = zeros(num_clusters,2);

% loop through clusters
for i = 1:num_clusters
   
   iclust = find(bins==i);   % cells in cluster i
   for iclust_i = iclust
      atag_mat(ctag_mat==iclust_i) = i; 
   end
   
   % calculate traction forces
   [ycluster, xcluster] = find(atag_mat==i); % all voxels in cluster i
   
   % cluster centroid
   y_cent = mean(ycluster);
   x_cent = mean(xcluster);

   cluster_centroids(i,:) = [x_cent y_cent];
   
   for k = 1:length(iclust)
      [y_cell,x_cell] = find(ctag_mat==iclust(k));
      fxnet(iclust(k)) = sum(x_cent-x_cell);  % net traction force for cell i
      fynet(iclust(k)) = sum(y_cent-y_cell);
      cell_centroids(iclust(k),:) = [mean(x_cell) mean(y_cell)]; % centroid of each cell
      
      % collect in a vector for quiver plot
      traction_x = [traction_x; x_cell];
      traction_y = [traction_y; y_cell];
      traction_u = [traction_u; x_cent-x_cell];
      traction_v = [traction_v; y_cent-y_cell];
   end
       
end

%JUNCTION MATRIX : RANK DEFICIENT
[m,n] = find(Aconn);
J = sparse(m, (m-1)*Ncell + n, 1, Ncell, Ncell^2);
% define coefficients for the following junction-traction force equations:
% sum J_ij = - fnet_i, where j are the indices of cells connected to cell i

%ASSUMPTION MATRIX OF REACTION FORCES
Njunc = Ncell*(Ncell-1)/2;
connu=triu(Aconn,1); %upper triangle of connectivity matrix
[i,j]=find(connu);      %nonzero indices of junctions
p=j+(i-1)*Ncell;          
q=i+(j-1)*Ncell;          
k=1:length(p);
A = sparse(k,p,1,Njunc, Ncell^2) + sparse(k,q,1,Njunc, Ncell^2);
% defines coefficients for the following junction force symmetry equations:
% J_ij + J_ji = 0

J=[J;A];
% T = [-fxnet -fynet; zeros(Njunc,2)];
% x = J\T;
% jx = x(:,1);  % Ncell^2 x 1 vector 
% jy = x(:,2);

T = [-fxnet; zeros(Njunc,1)];
% jx = J\T;
 jx = lsqminnorm(J,T);
% jx = (J'*J)\(J'*T);

T = [-fynet; zeros(Njunc,1)];
% jy = J\T;
 jy = lsqminnorm(J,T);
% jy = (J'*J)\(J'*T);

% Ncell x Ncell array
jx_mat = reshape(jx, Ncell, Ncell)'; % junction force between cell i and j
jy_mat = reshape(jy, Ncell, Ncell)';
neighbors_mat = zeros(Ncell);
jmag_mat = sqrt(jx_mat.^2 + jy_mat.^2);

% junc forces distributed across cell boundaries
[m,n] = find(Aconn);
jx_node_mat = zeros(ny, nx);
jy_node_mat = zeros(ny, nx);

% collect (x, y, jx, jy) for each voxel (otherwise 0)
% "norm" junction forces are normalized by # of boundary points or
% neighbors
xjunc_vox = [];
yjunc_vox = [];
jx_vox = [];
jy_vox = [];
jx_vox_norm = [];
jy_vox_norm = [];
for i = 1:length(m)
   [Nb, xb, yb] = find_neighbors(ctag_mat, n(i), m(i));
   neighbors_mat(n(i), m(i)) = Nb;
   for q = 1:Nb
      jx_node_mat(yb(q), xb(q)) =  jx_mat(n(i), m(i))/Nb;
      jy_node_mat(yb(q), xb(q)) =  jy_mat(n(i), m(i))/Nb;
      
      xjunc_vox = [xjunc_vox; xb(q)];
      yjunc_vox = [yjunc_vox; yb(q)];
      jx_vox_norm = [jx_vox_norm; jx_mat(n(i), m(i))/Nb];
      jy_vox_norm = [jy_vox_norm; jy_mat(n(i), m(i))/Nb];
      jx_vox = [jx_vox; jx_mat(n(i), m(i))];
      jy_vox = [jy_vox; jy_mat(n(i), m(i))];
   end
end

%  junc force magnitude at midpoint b/w cell centroids
[m,n] = find(Aconn);

% collect (x, y, jx, jy), for each junction where (x,y) is the center
% between cell centroids
xjunc_cent = [];
yjunc_cent = [];
jx_cent = [];
jy_cent = [];
jmag_cent = [];
for i = 1:length(m)
   [y1,x1] = find(ctag_mat==m(i));
   [y2,x2] = find(ctag_mat==n(i));
   
   yjunc = (mean(y1)+mean(y2))/2;
   xjunc = (mean(x1)+mean(x2))/2;
   
   xjunc_cent = [xjunc_cent; xjunc];
   yjunc_cent = [yjunc_cent; yjunc];
   jx_cent = [jx_cent; jx_mat(n(i), m(i))];
   jy_cent = [jy_cent; jy_mat(n(i), m(i))];
   jmag_cent = [jmag_cent; jmag_mat(n(i), m(i))];
end

   

if 0
   imagesc(atag_mat); hold on;
   plot(centroids(:,1), centroids(:,2),'ko','linewidth',2);
   quiver(traction_x, traction_y, traction_u, traction_v, 'w');
end
   
if 0
    figure;
   imagesc(ctag_mat); hold on;
   plot(centroids(:,1), centroids(:,2),'ko','linewidth',2);
   quiver(cell_centroids(:,1), cell_centroids(:,2),...
       fxnet, fynet, 'w');
    quiver(x_node_all, y_node_all, jx_node_all, jy_node_all,'k');
%     quiver(xjunc_all, yjunc_all, jx_all, jy_all, 'k');

end

if 0
    figure
    imagesc(sqrt(jx_node_mat.^2 + jy_node_mat.^2));
end


1;
end

function [Nb, xb, yb] = find_neighbors(ctag_mat, m, n)
% code to find the voxels in cell m that "neighbor" with cell n 
% by design, the connection matrix Aconn(m,n) = Aconn(n,m) = 1 
% (otherwise no border b/w these cells)
% Note one voxel can be included multiple in the outputs, if the cells share
% a border on multiple sides of the voxel

[y1, x1] = find(ctag_mat == m);
z = zeros(size(ctag_mat));
z(ctag_mat == m) = 1;
z(ctag_mat == n) = 2;

xb = []; yb = [];

for p = 1:length(x1)
    xi = x1(p); yi = y1(p);  
    % loop through all 4 neighbors
    for i = 1:4
       switch i
           case 1
               xn = xi-1; yn = yi;
           case 2
               xn = xi+1; yn = yi;
           case 3
               xn = xi; yn = yi-1;
           case 4
               xn = xi; yn = yi+1;
       end
       if z(yn,xn) == 2
           xb = [xb xi];
           yb = [yb yi];
       end
    end
    
    
    
end
Nb = length(xb);
end
   
