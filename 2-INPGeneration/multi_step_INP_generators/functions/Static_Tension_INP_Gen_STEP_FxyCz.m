%% Fiber Network Simulation, INP file generation main script
% function [] = Static_Tension_INP_Gen(file_loc,save_loc,strain)

function [final_nodes, final_fibers] = Static_Tension_INP_Gen_STEP_FxyCz(file_loc,strain, l_fiber)   
load_file=file_loc;
load(load_file,'final_nodes','final_fibers')
% 2ef
[final_fibers, final_nodes] = add_midnode_in_fiber_3D(final_fibers, final_nodes) ;

side_lengthx=round((max(final_nodes(:,2))-min(final_nodes(:,2))),-1);
side_lengthy=round((max(final_nodes(:,3))-min(final_nodes(:,3))),-1);
side_length_z = max(final_nodes(:,4))-min(final_nodes(:,4)) ;
% delta = side_length*strain/2;
delta = side_lengthy*strain;
delta_z = side_length_z * strain ;
strain_str=num2str(strain);
% file = save_loc;
% file = [save_loc,'.inp'];
file = [file_loc,'_strain_',strain_str,'.inp'];


% MATERIAL PROPERTIES
% Ratio of bending stiffness to stretching stiffness:
% Stiffness ratio = SR = EI / EA*l_fiber^2
SR = 10^(-4);  % THIS IS THE MODEL'S DIMENSIONLESS STIFFNESS %%10-4 standard%%
E = 1; % Set Young's modulus to 1 so stress is nondimensionalized. Note 
% that area is also set to 1, so the quantity EA equals 1. Thus all
% concentrated forces (cload) are in units of EA.
nu = 0;   % Poisson's ratio. If we don't give Abaqus a value, it uses 0.
% Fiber cross sectional area
Area = 1;
% l_fiber = 0.0002 ; % added by Mainak
% Fiber moment of inertia. Note that SR = EI / EA*lc^2
MOI = SR*Area*(l_fiber)^2;
% MOI = (Area^2)/(4*pi);
J=2*MOI;

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


%%============== Finding Edges ============================================
edge_x_pos = find(final_nodes(:,2) > (side_lengthx/2 - 25) );
edge_x_neg = find(final_nodes(:,2) < (-1*side_lengthx/2 + 25) );
edge_y_pos = find(final_nodes(:,3) > (side_lengthy/2 - 120) );
edge_y_neg = find(final_nodes(:,3) < (-1*side_lengthy/2 + 30) );

edge_z_pos = find(final_nodes(:,4) > (side_length_z/2 - 2.25) );
edge_z_neg = find(final_nodes(:,4) < (-1*side_length_z/2 + 2.25 ) );
side_length_z ;


    
%%
midpoint_nodes = zeros(length(final_fibers),4);

for j=1:length(final_fibers)
    midpt = .5*(final_nodes(final_fibers(j,2),2:4)+final_nodes(final_fibers(j,3),2:4));
    midpoint_nodes(j,:) = [length(final_nodes)+j midpt];
end

% Redo the elements to include the midpoint nodes
final_fibers = [final_fibers(:,1) final_fibers(:,2) midpoint_nodes(:,1) final_fibers(:,3)];

% Add midpoint nodes to final_nodes
final_nodes = [final_nodes; midpoint_nodes];


%% --- WRITE OUTPUT TO DATA FILE ---

fid = fopen(file,'w');

fprintf(fid,'*HEADING\n');
fprintf(fid,'**\n');
fprintf(fid,'** 3D TENSION Model Definition\n');
fprintf(fid,'**\n');


% --- NODES ---
fprintf(fid,'*NODE, NSET=Nodes\n');
fprintf(fid,'%8.0f,%20.12E,%20.12E,%20.12E\n',final_nodes');
fprintf(fid,'*NSET, NSET=XPOS\n');
fprintf(fid,'%8.0f\n',edge_x_pos');
fprintf(fid,'*NSET, NSET=XNEG\n');
fprintf(fid,'%8.0f\n',edge_x_neg');
fprintf(fid,'*NSET, NSET=YPOS\n');
fprintf(fid,'%8.0f\n',edge_y_pos');
fprintf(fid,'*NSET, NSET=YNEG\n');
fprintf(fid,'%8.0f\n',edge_y_neg');
fprintf(fid,'*NSET, NSET=ZPOS\n');
fprintf(fid,'%8.0f\n',edge_z_pos');
fprintf(fid,'*NSET, NSET=ZNEG\n');
fprintf(fid,'%8.0f\n',edge_z_neg');


% --- ELEMENTS ---
fprintf(fid,'**\n');
fprintf(fid,'** Elements\n');
fprintf(fid,'**\n');

%linear beam is B31 in 3D, B21 in 2D
%Quadratic beam is B32 in 3D, B22 in 2D
fprintf(fid,'*ELEMENT, TYPE=B32, ELSET=EL_fiber\n'); % For beam elements
fprintf(fid,'%8.0f,%8.0f,%8.0f,%8.0f\n',final_fibers');


% % --- SECTION ASSIGNMENT AND MATERIAL PROPS ---

fprintf(fid,'**\n');
fprintf(fid,'** Section Assignment and Material Properties\n');
fprintf(fid,'**\n');

fprintf(fid,'*BEAM GENERAL SECTION, ELSET=EL_fiber, DENSITY=%20.10E, SECTION=GENERAL\n',density); 
% No material given (it's defined here). Density is a required parameter.
% Poisson's ratio is optional; leaving blank sets it to zero.
fprintf(fid,'%20.10E,%20.10E,%20.10E,%20.10E,%20.10E\n',Area,MOI,0,MOI,J); % Area | I11 | I12 | I22 | J
fprintf(fid,'\n'); % blank line to accept default direction consines
fprintf(fid,'%20.10E,%20.10E\n',E,E/(2*(1+nu))); % Young's modulus | shear modulus | other entries unneeded
fprintf(fid,'*SECTION POINTS\n'); % Section points are needed to compute stress and strain within element
fprintf(fid,'%10.4f,%10.4f\n',0.5,0); % This gives x1 and x2 position of section points
fprintf(fid,'*DAMPING, ALPHA=%20.10E\n',alpha_r);

%---Load Steps---
fprintf(fid,'**\n');
fprintf(fid,'** Load Steps\n');
fprintf(fid,'**\n');

delta = delta / 20 ; % Twenty equal load steps

% --- STEP 1-20 ---
for ik = 1:20
    
fprintf(fid,'*STEP,NAME=STEP00%.3d, NLGEOM=YES, INC=300000\n',ik);
fprintf(fid,'*DYNAMIC, APPLICATION=QUASI-STATIC\n');
fprintf(fid,'%10.4f,%8.0f,%20.3e,%8.2f\n',[dt duration 0 0]); 


%     fprintf(fid,'*STEP,NAME=STEP00%.3d, NLGEOM=YES, INC=1000\n',1);
%     fprintf(fid,'*STATIC\n');
%     fprintf(fid,'%10.4f,%8.0f,%20.3e,%8.2f\n',[dt duration 0 0]); 


    fprintf(fid,'*BOUNDARY, TYPE=DISPLACEMENT\n');
    fprintf(fid,'YPOS, 2, , %8.4f\n', delta * ik);
    fprintf(fid,'*BOUNDARY, TYPE=DISPLACEMENT\n');
    fprintf(fid,'YNEG, 2, , %8.4f\n', 0);
%     fprintf(fid,'*BOUNDARY, TYPE=DISPLACEMENT\n');
%     fprintf(fid,'YNEG, 1, , %8.4f\n', 0);
    
    fprintf(fid,'*BOUNDARY, TYPE=DISPLACEMENT\n');
    fprintf(fid,'ZPOS, 3, , %8.4f\n', 0);
%     fprintf(fid,'*BOUNDARY, TYPE=DISPLACEMENT\n');
%     fprintf(fid,'XPOS, 2, , %8.4f\n', 0);
    fprintf(fid,'*BOUNDARY, TYPE=DISPLACEMENT\n');
    fprintf(fid,'ZNEG, 3, , %8.4f\n', 0);
%     fprintf(fid,'*BOUNDARY, TYPE=DISPLACEMENT\n');
%     fprintf(fid,'XNEG, 2, , %8.4f\n', 0);
     
%     fprintf(fid,'*BOUNDARY, TYPE=DISPLACEMENT\n');
%     fprintf(fid,'YPOS, 1, , %8.4f\n', delta/2);
%     fprintf(fid,'*BOUNDARY, TYPE=DISPLACEMENT\n');
% %   fprintf(fid,'YNEG, 2, , %8.4f\n', -delta);
%     fprintf(fid,'YNEG, 1, , %8.4f\n', -delta/2);

    
%     fprintf(fid,'*BOUNDARY, TYPE=DISPLACEMENT\n');
%     fprintf(fid,'ZPOS, 3, , %8.4f\n', delta_z);
%     fprintf(fid,'*BOUNDARY, TYPE=DISPLACEMENT\n');
% %   fprintf(fid,'YNEG, 2, , %8.4f\n', -delta);
%     fprintf(fid,'ZNEG, 3, , %8.4f\n', 0);
    

    fprintf(fid,'*OUTPUT, FIELD, NUMBER INTERVAL=20 \n');
    fprintf(fid,'*NODE OUTPUT\n');
    fprintf(fid,'U, RF\n');
    fprintf(fid,'*ELEMENT OUTPUT\n');
    fprintf(fid,'SE, LE, S, SK, SF, NFORC\n'); % SK gives section curvature change
    fprintf(fid,'*OUTPUT, HISTORY, VARIABLE=ALL, NUMBER INTERVAL=20\n');

    fprintf(fid,'*END STEP\n');

end
fclose('all');

return
