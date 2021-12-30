%% Fiber Network Simulation, INP file generation main script

% This script takes the workspace output from
function [] = Static_2D_Cyl_Ind_INP(file_loc, w, R_IND)

% Fiber_Network_Creation_Main_Script and generates an Abaqus .inp file for
% the FEA model

load_file=[file_loc,'.mat'];
load(load_file) %#ok<LOAD>

Lx=round((max(final_nodes(:,2))-min(final_nodes(:,2)))); %#ok<*NODEF>
Ly=round((max(final_nodes(:,3))-min(final_nodes(:,3))));


file = [file_loc,'_CylPunch_R_',num2str(R_IND),'.inp'];

% MATERIAL PROPERTIES

% Ratio of bending stiffness to stretching stiffness:
% Stiffness ratio = SR = EI / EA*l_fiber^2
SR = 10^(-4);  % THIS IS THE MODEL'S DIMENSIONLESS STIFFNESS
E = 1; % Set Young's modulus to 1 so stress is nondimensionalized. Note 
% that area is also set to 1, so the quantity EA equals 1. Thus all
% concentrated forces (cload) are in units of EA.
nu = 0;   % Poisson's ratio. If we don't give Abaqus a value, it uses 0.
% Fiber cross sectional area
Area = 1;
% Fiber moment of inertia. Note that SR = EI / EA*lc^2
MOI = SR*Area*l_fiber^2;
% MOI = (Area^2)/(4*pi);
J=2*MOI;

% Dynamic Parameters:
dt=8e-4;    % time increment
% Analysis duration
% Time over which to ramp the force
ramp_time_percent = .5;
% % Bulk viscosity. Used to prevent ringing.
% viscosity = 1e-2; % Default is 0.06
alpha_r = .007;  % damping
% Density for dynamic analysis
density = .1;


% this is for if you have multiple steps in your job
%duration_vector = [2500 2500 2500 4000];% 1500 1600 1800 2000];
%duration = sum(duration_vector);
duration = 2000;
ramp_time = ramp_time_percent*duration;

% everything is defined in the global cartesian CS, so for a rotation of
% the top surface we need to solve for what the cartesian displacements of
% each of the cload_nodes should be.

% Need to find the theta unit vector for each cload_node expressed in x and
% y:


%%============== Finding Edges ============================================
edge_x_pos = find(final_nodes(:,2) > (Lx/2 - 0.5) );
edge_x_neg = find(final_nodes(:,2) < (-1*Lx/2 + 0.5) );
edge_y_pos = find(final_nodes(:,3) > (Ly/2 - 0.5) );
edge_y_neg = find(final_nodes(:,3) < (-1*Ly/2 + 0.5) );
Contact_Nodes = final_nodes(boundary(final_nodes(:,2),final_nodes(:,3)),1);
    
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
final_nodes(:,4)=0;


%% Create Indentor 
R=R_IND;
xo=0; yo=Ly/2+R;
xa=R; ya=yo;
xb=0; yb=yo-R;
xc=-R; yc=yo;
RI_id = final_nodes(end,1)+1;
indent_nodes = [RI_id, xo,yo,0;...
                RI_id+1, xa,ya,0;...
                RI_id+2, xb,yb,0;...
                RI_id+3, xc,yc,0];

%% --- WRITE OUTPUT TO DATA FILE ---

fid = fopen(file,'w');

fprintf(fid,'*HEADING\n');
fprintf(fid,'**\n');
fprintf(fid,'** 2D TENSION Model Definition\n');
fprintf(fid,'**\n');
fprintf(fid,'**Info about this job:\n');
fprintf(fid,'**********************************\n');
fprintf(fid,'**SR = %8.1E\n',SR);
fprintf(fid,'**domain diameter = %8.1E\n',lx);
fprintf(fid,'**N_fibers = %8.1E\n',length(final_fibers));
fprintf(fid,'**alpha (damping) = %8.1E\n',alpha_r);
% fprintf(fid,'**rotation angle = %3.2f\n',torsion_angle);
fprintf(fid,'**total duration = %8.2E\n', duration);
fprintf(fid,'**ramp time fraction = %8.2E\n',ramp_time_percent);
fprintf(fid,'**timestep = %8.1E\n',dt);

fprintf(fid,'*PREPRINT,HISTORY=NO,MODEL=NO\n');
fprintf(fid,'**\n');


% --- NODES ---
fprintf(fid,'*NODE, NSET=Nodes\n');
fprintf(fid,'%8.0f,%20.12E,%20.12E,%20.12E\n',final_nodes');
fprintf(fid,'*NODE, NSET=Indentor\n');
fprintf(fid,'%8.0f,%20.12E,%20.12E,%20.12E\n',indent_nodes');
fprintf(fid,'*NSET, NSET=XPOS\n');
fprintf(fid,'%8.0f\n',edge_x_pos');
fprintf(fid,'*NSET, NSET=XNEG\n');
fprintf(fid,'%8.0f\n',edge_x_neg');
fprintf(fid,'*NSET, NSET=YPOS\n');
fprintf(fid,'%8.0f\n',edge_y_pos');
fprintf(fid,'*NSET, NSET=YNEG\n');
fprintf(fid,'%8.0f\n',edge_y_neg');
fprintf(fid,'*NSET, NSET=CONTACT\n');
fprintf(fid,'%8.0f\n',Contact_Nodes');


% --- ELEMENTS ---
fprintf(fid,'**\n');
fprintf(fid,'** Elements\n');
fprintf(fid,'**\n');
%linear beam is B31 in 3D, B21 in 2D
%Quadratic beam is B32 in 3D, B22 in 2D
fprintf(fid,'*ELEMENT, TYPE=B22, ELSET=EL_fiber\n'); % For beam elements
fprintf(fid,'%8.0f,%8.0f,%8.0f,%8.0f\n',final_fibers');

% --- SURFACES ---
fprintf(fid,'*SURFACE, TYPE=SEGMENTS, NAME=INDENT\n');
    fprintf(fid,'START, %20.10E,%20.10E\n',xa,ya);
    fprintf(fid,'CIRCL, %20.10E,%20.10E,%20.10E,%20.10E\n',xb,yb,xo,yo);
    fprintf(fid,'CIRCL, %20.10E,%20.10E,%20.10E,%20.10E\n',xc,yc,xo,yo);
    fprintf(fid,'*RIGID BODY, ANALYTICAL SURFACE=INDENT, REF NODE=%8.0f\n',RI_id);
fprintf(fid,'*SURFACE, TYPE=NODE, NAME=CONTACT\n');
    fprintf(fid,'CONTACT\n');
    


% CREATE CONTACT PAIR
fprintf(fid,'*Surface Interaction,NAME=CONTACT1\n');
fprintf(fid,'*CONTACT PAIR,INTERACTION=CONTACT1\n');
fprintf(fid,'CONTACT,INDENT\n');



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
    

%--- BOUNDARY CONDITIONS ---
fprintf(fid,'**\n');
fprintf(fid,'** Boundary Conditions\n');
fprintf(fid,'**\n');

% Fixed
fprintf(fid,'*BOUNDARY\n');
% fprintf(fid,'XNEG, ENCASTRE\n');
% fprintf(fid,'XPOS, ENCASTRE\n');
fprintf(fid,'YNEG, ENCASTRE\n');

%---Load Steps---
fprintf(fid,'**\n');
fprintf(fid,'** Load Steps\n');
fprintf(fid,'**\n');


fprintf(fid,'*STEP,NAME=STEP00%.3d, NLGEOM=YES, INC=1000\n',1);
fprintf(fid,'*STATIC\n');
fprintf(fid,'%10.4f,%8.0f,%20.3e,%8.2f\n',[dt duration 0 0]); 
% suggested initial time increment | time period | min time increment | max time icrement
% zeros allow for abaqus defaults



fprintf(fid,'*BOUNDARY, TYPE=DISPLACEMENT\n');
fprintf(fid,'Indentor, 1, , %8.4f\n', 0);
fprintf(fid,'*BOUNDARY, TYPE=DISPLACEMENT\n');
fprintf(fid,'Indentor, 2, , %8.4f\n', -w);
fprintf(fid,'*BOUNDARY, TYPE=DISPLACEMENT\n');
fprintf(fid,'Indentor, 4,6, %8.4f\n', 0);


fprintf(fid,'*OUTPUT, FIELD, NUMBER INTERVAL=20 \n');
fprintf(fid,'*NODE OUTPUT\n');
fprintf(fid,'U, V, A, RF\n');
fprintf(fid,'*ELEMENT OUTPUT\n');
fprintf(fid,'SE, LE, S, SK, SF, NFORC\n'); % SK gives section curvature change
fprintf(fid,'*CONTACT OUTPUT, SURFACE=CONTACT, SECOND SURFACE=INDENT\n');
fprintf(fid,'CNAREA, CSTRESS, CDISP, CFORCE, CNAREA\n');
fprintf(fid,'*OUTPUT, HISTORY, VARIABLE=ALL, NUMBER INTERVAL=20\n');
fprintf(fid,'*END STEP\n');

        
fclose('all');