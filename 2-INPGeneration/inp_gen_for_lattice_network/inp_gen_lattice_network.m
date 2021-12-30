% define material and generate an Abaqus input file for lattice fiber
% network
% Several sections are modified in the existing INP generator of our group by Mainak Sarkar


function [final_fibers, final_nodes] = inp_gen_lattice_network(file, strain, l_fiber, final_fibers, final_nodes, SR, edge_thkness)

% make the folder that will contain the INP file of network
% P = 'C:\Users\msarkar3\Desktop\MATLAB_Path\3D_printing_2D';
% F = 'lattice_fiber';
% mkdir(P,F)
% new = fullfile(P,F) ;

side_length=round((max(final_nodes(:,2))-min(final_nodes(:,2))),-1) ;
delta = side_length*strain/2 ;
% strain_str=num2str(strain);
% file = save_loc;
% file = [save_loc,'.inp'];
% file = [new,'\strain_',strain_str,'.inp'] ;

% MATERIAL PROPERTIES
% Ratio of bending stiffness to stretching stiffness:
% Stiffness ratio = SR = EI / EA*l_fiber^2
% SR = 10^(-4);  % THIS IS THE MODEL'S DIMENSIONLESS STIFFNESS %%10-4 standard%%
E = 1; % Set Young's modulus to 1 so stress is nondimensionalized. Note 
% that area is also set to 1, so the quantity EA equals 1. Thus all
% concentrated forces (cload) are in units of EA.
nu = 0;   % Poisson's ratio. If we don't give Abaqus a value, it uses 0.
% Fiber cross sectional area
Area_a = 1;
% l_fiber = 0.0002 ; % added by Mainak
% Fiber moment of inertia. Note that SR = EI / EA*lc^2
MOI_a = SR*Area_a*(l_fiber)^2;
% MOI = (Area^2)/(4*pi);
J_a = 2*MOI_a ;


% Parameters:
dt=8e-4;    % time increment
alpha_r = .007;  % damping
% Density for dynamic analysis
density = .1;
duration = 2000;


% everything is defined in the global cartesian CS, so for a rotation of
% the top surface we need to solve for what the cartesian displacements of
% each of the cload_nodes should be.

% Need to find the theta unit vector for each cload_node expressed in x and
% y:

side_length = max(final_nodes(:,3)) - min(final_nodes(:,3)) ; 

d = edge_thkness * side_length ;
%%============== Finding Edges (do edits if necessary)============================================
edge_x_pos = find(final_nodes(:,2) > (side_length - 8) ) ;
edge_x_neg = find(final_nodes(:,2) < ( + 8) ) ;
edge_y_pos = find(final_nodes(:,3) > (side_length - d) ) ;
edge_y_neg = find(final_nodes(:,3) < ( + d) ) ;

    
%%
midpoint_nodes = zeros(length(final_fibers),3);

for j=1:length(final_fibers)
    midpt = .5*(final_nodes(final_fibers(j,2),2:3)+final_nodes(final_fibers(j,3),2:3));
    midpoint_nodes(j,:) = [length(final_nodes)+j midpt];
end

% Redo the elements to include the midpoint nodes
final_fibers = [final_fibers(:,1) final_fibers(:,2) midpoint_nodes(:,1) final_fibers(:,3)] ;
% [final_fibers_a, final_fibers_b, final_fibers_c, final_fibers_d] = group_maker(fibers, nodes, r1, r2, r3, final_fibers, O1, O2) ;
% Add midpoint nodes to final_nodes
final_nodes = [final_nodes; midpoint_nodes];


%% --- WRITE OUTPUT TO DATA FILE ---

fid = fopen(file,'w');

fprintf(fid,'*HEADING\n');
fprintf(fid,'**\n');
fprintf(fid,'** 2D TENSION Model Definition\n');
fprintf(fid,'**\n');


% --- NODES ---
fprintf(fid,'*NODE, NSET=Nodes\n');
% fprintf(fid,'%8.0f,%20.12E,%20.12E,%20.12E\n',final_nodes');
fprintf(fid,'%8.0f,%20.12E,%20.12E\n',final_nodes');

fprintf(fid,'*NSET, NSET=XPOS\n');
fprintf(fid,'%8.0f\n',edge_x_pos');
fprintf(fid,'*NSET, NSET=XNEG\n');
fprintf(fid,'%8.0f\n',edge_x_neg');
fprintf(fid,'*NSET, NSET=YPOS\n');
fprintf(fid,'%8.0f\n',edge_y_pos');
fprintf(fid,'*NSET, NSET=YNEG\n');
fprintf(fid,'%8.0f\n',edge_y_neg');


% --- ELEMENTS ---
fprintf(fid,'**\n');
fprintf(fid,'** Elements\n');
fprintf(fid,'**\n');

%linear beam is B31 in 3D, B21 in 2D
%Quadratic beam is B32 in 3D, B22 in 2D
fprintf(fid,'*ELEMENT, TYPE=B32, ELSET=EL_fiber_a\n'); % For beam element set a
% fprintf(fid,'%8.0f,%8.0f,%8.0f,%8.0f\n',final_fibers');
fprintf(fid,'%8.0f,%8.0f,%8.0f,%8.0f\n',final_fibers');

% % --- SECTION ASSIGNMENT AND MATERIAL PROPS ---

fprintf(fid,'**\n');
fprintf(fid,'** Section Assignment and Material Properties\n');
fprintf(fid,'**\n');

% section for element set a
fprintf(fid,'*BEAM GENERAL SECTION, ELSET=EL_fiber_a, DENSITY=%20.10E, SECTION=GENERAL\n',density); 
% No material given (it's defined here). Density is a required parameter.
% Poisson's ratio is optional; leaving blank sets it to zero.
fprintf(fid,'%20.10E,%20.10E,%20.10E,%20.10E,%20.10E\n',Area_a,MOI_a,0,MOI_a,J_a); % Area | I11 | I12 | I22 | J
fprintf(fid,'\n'); % blank line to accept default direction consines
fprintf(fid,'%20.10E,%20.10E\n',E,E/(2*(1+nu))); % Young's modulus | shear modulus | other entries unneeded
fprintf(fid,'*SECTION POINTS\n'); % Section points are needed to compute stress and strain within element
fprintf(fid,'%10.4f,%10.4f\n',0.5,0); % This gives x1 and x2 position of section points
fprintf(fid,'*DAMPING, ALPHA=%20.10E\n',alpha_r);
% end of section assignment for element set a


%---Load Steps---
fprintf(fid,'**\n');
fprintf(fid,'** Load Steps\n');
fprintf(fid,'**\n');


% --- STEP 1 ---
fprintf(fid,'*STEP,NAME=STEP00%.3d, NLGEOM=YES, INC=300\n',1);
fprintf(fid,'*DYNAMIC, APPLICATION=QUASI-STATIC\n');
fprintf(fid,'%10.4f,%8.0f,%20.3e,%8.2f\n',[dt duration 0 0]); 

%{
    fprintf(fid,'*STEP,NAME=STEP00%.3d, NLGEOM=YES, INC=1000\n',1);
    fprintf(fid,'*STATIC\n');
    fprintf(fid,'%10.4f,%8.0f,%20.3e,%8.2f\n',[dt duration 0 0]); 
%}

    fprintf(fid,'*BOUNDARY, TYPE=DISPLACEMENT\n');
    fprintf(fid,'YPOS, 2, , %8.4f\n', delta);
    fprintf(fid,'*BOUNDARY, TYPE=DISPLACEMENT\n');
    fprintf(fid,'YNEG, 2, , %8.4f\n', -delta);

    fprintf(fid,'*OUTPUT, FIELD, NUMBER INTERVAL=20 \n');
    fprintf(fid,'*NODE OUTPUT\n');
    fprintf(fid,'U, RF\n');
    fprintf(fid,'*ELEMENT OUTPUT\n');
    fprintf(fid,'SE, LE, S, SK, SF, NFORC\n'); % SK gives section curvature change
    fprintf(fid,'*OUTPUT, HISTORY, VARIABLE=ALL, NUMBER INTERVAL=20\n');

    fprintf(fid,'*END STEP\n');

        
fclose('all');