% 2D lattice-based fiber network generator code
% Consists of 8 sections.
% Section 1 is taken from Prof Jacob Notbohm's code with some modifications
% Sections 2 to 8 in this main file and all the supporting function files
% are written by Mainak Sarkar, University of Wisconsin-Madison
tic
clear all
clc

%% Define inputs: INPUT PANEL
% Domain dimensions
lx = 110 ; 
ly = 110 ;
rho_nodal = 0.036 ; % input desired nodal density of network / For triangulation, method 2 (pre-midnode addn.).
% But, for VT, method 2, actual density (pre-midnode) = rho_nodal * (~1.8)
% probability of existance of a fiber between two neighboring nodes
pr = 0.97 ; % tested values 0.92, 0.925, ..., 1

% enter the anticipated rectangular cross-section of fiber: x-y width = w
% and z-height = h
w = 0.5 ; 
h = 1.5 ;
% specify required mesh
% Method 1: Tri-lattice base with each node connected to exacly 6 fibers.
% 
% A special note: A hexagonal mesh is extracted from this (optional) to generate a network
% consisting of hexagonal lattice elements with the condition that each
% node is connected to exactly 3 fibers.

% Method 2: Tri-lattice developed with delaunay triangulation considering
% our random nodes, but there is no constraint on the number of fibers that
% can be attached to each node to satisfy the concept of delaunay
% triangulation.

method_lattice_triangle = 2 ; % input 1 for method 1 and input 2 for method 2
dtvo = 1 ; % 0 for delaunay triangulation only**, 1 for DT plus voronoi diagram from DT (WORKS IF METHOD 2 IS ON)
require_hexagonal_mesh = 1 ; % 0 or 1, INACTIVE / OFF BY DEFAULT if method_2 is ON
need_edge_treatment = 1 ; % 0 or 1. If 1, then fiber remover will not operate near edges (specified by 'edge_thkness').
edge_thkness = 15/ly ; % what fraction of side length is the intended edge thickness (top and bottom edges)? 
% want a robust edge for polygonal lattice network?
robust_edge = 1 ; % 0 or 1. Choose 1 for polygonal lattice-based network only.

% make the folder that will contain the MAT file of network
date = datestr(now,'mm_dd_yyyy') ;
P = 'D:\';
F = ['random_fiber_net_2D_domain_',num2str(lx),'x',num2str(ly),'_rho_',num2str(rho_nodal),'_pr_',num2str(pr),'_w_',num2str(w),'_h_',num2str(h),'_date_', date] ; % in the MATLAB current folder
mkdir(P,F)
new = fullfile(P,F) ;

% Meshing parameters
[Nx, Ny, dNx, dNy] = find_node_numbers_2Dlattice(rho_nodal, lx, ly) ; % finds number of nodes in x and y dimension on domain

knode = 1 ;
for m=1:Ny % rows
    for n=1:Nx % columns
        % Node number
        N1(knode,1) = knode;
        N1(knode,4) = m; % row number
        N1(knode,5) = n; % column number. Row and column numbers are used in later logic.
        
        if rem(m,2)~=0 % For odd rows
            % X coordinate
            N1(knode,2) = dNx*(n-1);
            % Y coordinate
            N1(knode,3) = dNy*(m-1);
        else % For even rows
            % X coordinate
            N1(knode,2) = dNx*(n-1) + dNx/2; % shift x coords to right by a half node spacing
            % Y coordinate
            N1(knode,3) = dNy*(m-1);
        end
        %{
        % Check node location and if it's along circular boundary, include
        % it in a new nodeset
        tol = 0.85*mean([dNx dNy]); % Tolerance here is roughly half node-to-node spacing
        for k = 1:length(a)
            if abs( sqrt((N1(knode,2)-xc(k))^2 + (N1(knode,3)-yc(k))^2) - a(k)) < tol
                nodes_circ_k = Ncirc{k};
                nodes_circ_k = [nodes_circ_k ; N1(knode)];
                % Ncirc is a 1D cell with the k-th element equal to a 
                % vector of nodes on the boundary of the k-th circle
                Ncirc{k}=nodes_circ_k;
            end
        end
        %}
        
       
        if m~=1 && m~=Ny && n~=1 && n~=Nx % If not on edges, add randomness
            randx = (rand(1)-0.5)*dNx * 0.9;
            randy = (rand(1)-0.5)*dNy * 0.9;
            N1(knode,2) = N1(knode,2) + randx;
            N1(knode,3) = N1(knode,3) + randy;
        end
        
        knode = knode + 1;
    end
end
%--------------------------------------END of SECTION 1--------------------------------------------------------------------------------------------- 

%% This is section 2 describing METHOD- 1. This and all the subsequent sections, supporting functions are written by Mainak.
if method_lattice_triangle == 1
% For checking nodal coordinates
% figure
% plot(N1(:,2),N1(:,3),'o')
% axis equal

% find neighbors of each arbitrary node (not edge nodes)
p = 1 ; % p = 1 implies six fibers are connected to one node
i = 0 ; % loop counter 
for k = 2:1:(Ny-1)
for N = ((k-1)*Nx+2):1:(((k-1)*Nx+Nx)-1)
    i = i + 1 ;
[Neighbor_nodes] = arbit_node_neighbors(N, N1) ;
el_set = nodal_disorder(p, Neighbor_nodes, N1, N) ;
Cell_elements{i} = el_set ;
end
end
el_set_final = cell2mat(Cell_elements') ;
% size(el_set_final)

%{
size(el_set_final)
el_set_final = el_set_final(:,[2 1]) ;
el_set_final = unique(el_set_final,'rows') ;
el_set_final = el_set_final(:,[2 1]) 
size(el_set_final)
%}

% edge nodes: finding neighbors
[el_set_1, el_set_2, el_set_3, el_set_4, el_set_5, el_set_6, el_set_7, el_set_8] = edge_element_generator(Nx, Ny, N1) ;
el_set_edges = [el_set_1; el_set_2; el_set_3; el_set_4; el_set_5; el_set_6; el_set_7; el_set_8] ;

% conglomerate all fibers
el_set_final = [el_set_final ; el_set_edges] ;

%{
% eliminate repeated fibers
for kk = 1:size(el_set_final,1) 
    for ll = 1:size(el_set_final,1)
    if el_set_final(kk,1) == el_set_final(ll,2)
        if el_set_final(ll,1) == el_set_final(kk,2)
           el_set_final(kk,:) = zeros(1, size(el_set_final,2)) ;
        end
    end
    end
end
%}  

% el_set_final = el_set_final(any(el_set_final,2),:) ;
% size(el_set_final) 

% eliminate repeated fibers (more efficient)
el_set_final = unique(sort(el_set_final,2),'rows') ;


% fiber consideration on particular node/ No of fibers linked to the
% selected arbit node, N

el_set_final = [(1:size(el_set_final,1))' el_set_final] ;
end
nodes_set_final = N1(:,1:3) ;

%{
%% deleting fibers from each node (optional)
prob = 4/6 ;
if prob < 1 
[el_set_final, nodes_set_final] = fiber_delete_assist(prob, Nx, Ny, N1, el_set_final, nodes_set_final) ;
end
%}

%{
%% randomly removing fibers from the entire domain
[el_set_final] = random_fiber_eliminator(pr, el_set_final) ;  % use it
% after hex mesh gen if that is present
%}

%----------------------------------------END OF METHOD 1----------------------------------------------------------------------------------------------------------------------

%% This is sction 3. It is an alternative, fast triangular mesh generator (METHOD- 2). If this method is ON, 
% it will replace the result of METHOD-1 and keep only the results of METHOD- 2. 
% Please refrain from switching on hex_mesh_gen if this is ON
if method_lattice_triangle == 2
 [nodes_set_final, el_set_final] = delaunay_tri_mesh_gen(nodes_set_final, dtvo)   ;
end


% Second part of method 2 (subsection 3.1): Obtain voronoi diagram from the
% delaulay triangulation is included in the delaunay_tri_mesh_gen function.



%----------------------------------------END OF METHOD 2-------------------------------------------------------------------------------------------------------

%% This is sction 4. It makes hexagonal mesh (random) by 'delaunay triangulation' (optional) (ONLY use it if METHOD- 1 is ON, METHOD- 2 is OFF)
if require_hexagonal_mesh == 1
if method_lattice_triangle ~= 2
    [el_set_final] = hexa_mesh_gen(N1, Nx, Ny, el_set_final, nodes_set_final) ;
end
end


%% This is sction 5. It randomly removes fibers from the entire domain (*may leave edges*) (**works universally in either methods**)
[el_set_final, nodes_set_final] = renumbering_el_nodes(el_set_final, nodes_set_final) ; % renumbering is necessary here
if need_edge_treatment == 0
[el_set_final] = random_fiber_eliminator(pr, el_set_final) ;
else
[el_set_final] = random_fiber_eliminator_adv(pr, el_set_final, nodes_set_final, edge_thkness) ;
end

%% This is sction 6. It removes dangling fibers (you can use this filter more than once if necessary) (**works universally in either methods**)
 [el_set_final, nodes_set_final] = dangling_fiber_remover(el_set_final, nodes_set_final) ;
 [el_set_final, nodes_set_final] = dangling_fiber_remover(el_set_final, nodes_set_final) ;
 [el_set_final, nodes_set_final] = dangling_fiber_remover(el_set_final, nodes_set_final) ;
 
%% (Optional section) If a robust edge is required for a polygonal lattice-based network.
if robust_edge == 1
[el_set_final, nodes_set_final] = renumbering_el_nodes(el_set_final, nodes_set_final) ;
[el_set_final] = draw_edge_fibers(nodes_set_final, el_set_final, edge_thkness) ;
end 

%% This is an addendum section 6(a). It eliminates ghost nodes (excess nodes not connected to any fibers).
 [nodes_set_final, node_index_refresher] = ghost_node_remover(el_set_final, nodes_set_final) ;

 
%% This is sction 7. It corrects indexing of fiber and node matrices (**works universally in either methods**)
[el_set_final, nodes_set_final] = renumbering_el_nodes(el_set_final, nodes_set_final) ; % this is necessary
[average_nodal_connectivity] = find_avg_node_connectivity(el_set_final, nodes_set_final, edge_thkness)    % average nodal connectivity at the lattice vertices, leave edges
% (away from edges)

%% This is sction 8. It reports anticipated value of kappa based on the c/s of fiber provided and also reports fiber length and nodal density of the network.
[l_fiber] = distance_finder(nodes_set_final, el_set_final, edge_thkness)  % get average fiber length 
I = h * (w^3)/12 ;
A = w * h ;
SR = (I/A)/(l_fiber)^2 % dimensionless bending stiffness ratio, kappa

ly_e = ly * (1-(2*edge_thkness)) ;
[effective_node_nos, rho_nodal_act] = rhonodal_finder(nodes_set_final, el_set_final, edge_thkness, lx, ly_e, 0) ;
rho_nodal_act

%% save the fibers and nodes of the network in a MAT file.
Path_lat = char(strcat(new,'\lattice_data_NM_',num2str(lx),'x',num2str(ly),'_rho_',num2str(rho_nodal),'.mat')) ;
save(Path_lat, 'nodes_set_final' , 'el_set_final', 'edge_thkness') ; 

toc

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
