% **This main script is an INP generator in tension boundary condition for use with the
% lattice-based fiber network generators
% **Note: Very minute modification is necessary in this main script and the corresponding function
% file to address other BCs
% written by Mainak Sarkar

%load the MAT file containing node and element data:
load('lattice_fiber_2020_07_07_sample_4.mat', 'el_set_final', 'nodes_set_final') 

%% Define tensile strain here along X:
% define strain (bc specific)
strain = 0.06 ; % enter tensile strain here (x-directional) [WILL CHANGE]
d = 0.2 ; % mention the diameter of fibers here to compute the dimensionless bending stiffness ratio, SR in next section

% make the folder that will contain the INP file of network
P = 'C:\Users\msarkar3\Desktop\MATLAB_Path\3D_printing_2D'; % give the path of your intended folder containing the INP file
F = 'lattice_fiber_2020_07_07_sample_4_inp'; % This is the name of the folder containing inp file
mkdir(P,F)
new = fullfile(P,F) ;
strain_str=num2str(strain);
% Define file names of the INP file
file = [new,'\strain_',strain_str,'.inp'] ;

edge_thkness = 0.05 ;
%% This section generates an inp file for static tensile displacement along x (**works universally in either methods**)
[l_fiber] = distance_finder(nodes_set_final, el_set_final) ; % get average fiber length 
SR = d^2/(16*(l_fiber)^2) ; % dimensionless bending stiffness ratio
[final_fibers, final_nodes] = inp_gen_lattice_network(file, strain, l_fiber, el_set_final, nodes_set_final, SR, edge_thkness) ;

