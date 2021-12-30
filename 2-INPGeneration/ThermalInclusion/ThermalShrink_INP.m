%Fiber Network Simulation, INP file generation main script
function ThermalShrink_INP(File_Path,Rad_particle,Rad_inclusion,Seed,Ang,Gap,Strain)
%% Load in File
load(File_Path,'final_nodes','final_fibers')
final_nodes(:,4)=0;

%% Find second inclusion to pin
R_nodal = sqrt(final_nodes(:,2).^2+final_nodes(:,3).^2);
R_min = round(min(R_nodal));
R_max = round(max(R_nodal));
Inc_R = Gap; %Radial Position of Center of Inclusion
Inc_T = Ang; %Angular Position of Center of Inclusion
    Inc_Td = round(Ang*180/pi);
[Inc_X,Inc_Y] = pol2cart(Inc_T,Inc_R);
Inc_diffx = final_nodes(:,2)-Inc_X;
Inc_diffy = final_nodes(:,3)-Inc_Y;
Inc_Nodes = final_nodes(sqrt(Inc_diffx.^2+Inc_diffy.^2) < Rad_inclusion,1);

%% Find Boundary Nodes
Outer_Nodes = final_nodes(R_nodal >= (R_max),1);
Shrink_Nodes = final_nodes(R_nodal <= (R_min),1);
    Shrink_Nodes_Full = final_nodes(R_nodal <= (R_min),:);

    
%% Create Inclusion
ElSize = 1.5;  %Assume lf = 1
[particle_nodes,particle_tri] = CreateInclusion(R_min,ElSize);
particle_nodes(:,1) = particle_nodes(:,1)+length(final_nodes);
particle_tri(:,1) = particle_tri(:,1)+length(final_fibers);
particle_tri(:,2:4) = particle_tri(:,2:4)+length(final_nodes);


%% File Name
File_Name = ['Domain',num2str(R_max),'_Strain',num2str(Strain),'_Inc',num2str(Rad_inclusion),...
                '-R',num2str(round(Inc_R),'%3.3d'),'T',num2str(Inc_Td,'%3.3d'),'_Seed',num2str(Seed),'.inp'];
File_Dir = fileparts(File_Path);
file = fullfile(File_Dir,File_Name);


%% MATERIAL PROPERTIES
% Ratio of bending stiffness to stretching stiffness:
% Stiffness ratio = SR = EI / EA*l_fiber^2
SR = 1e-4;  % THIS IS THE MODEL'S DIMENSIONLESS STIFFNESS
E = 1; % Set Young's modulus to 1 so stress is nondimensionalized. Note 
% that area is also set to 1, so the quantity EA equals 1. Thus all
% concentrated forces (cload) are in units of EA.
nu = 0;   % Poisson's ratio. If we don't give Abaqus a value, it uses 0.
% Fiber cross sectional area
Area = 1;
% Fiber moment of inertia. Note that SR = EI / EA*lc^2
Lfiber=1;
MOI = SR*Area*Lfiber^2;
J=MOI;

%% Dynamic Parameters:
dt=.1;    % time increment
% Analysis duration
duration = 2000;  % [s]
% Time over which to ramp the force
ramp_time = .75*duration;
% Density for dynamic analysis
density = .01;



%% --- WRITE OUTPUT TO DATA FILE ---

fid = fopen(file,'w');
fprintf(fid,'*HEADING\n');
fprintf(fid,'**\n');
fprintf(fid,'** Model Definition\n');
fprintf(fid,'**\n');
fprintf(fid,'**SR = %8.3E\n',SR);
fprintf(fid,'*PREPRINT,HISTORY=NO,MODEL=NO\n');
fprintf(fid,'**\n');


% --- NODES ---
fprintf(fid,'*NODE, NSET=Nodes\n');
fprintf(fid,'%8.0f,%20.12E,%20.12E,%20.12E\n',[final_nodes;particle_nodes]');

% --- NODE SETS ---
fprintf(fid,'*NSET, NSET=NBNDRY\n');
fprintf(fid,'%8.0f\n',Outer_Nodes');
fprintf(fid,'*NSET, NSET=CONNECT\n');
fprintf(fid,'%8.0f\n',Shrink_Nodes');
fprintf(fid,'*NSET, NSET=INCLUSION\n');
fprintf(fid,'%8.0f\n',Inc_Nodes');
fprintf(fid,'*NSET, NSET=PARTICLE\n');
fprintf(fid,'%8.0f\n',particle_nodes(:,1)');


% --- ELEMENTS ---
fprintf(fid,'**\n');
fprintf(fid,'** Elements\n');
fprintf(fid,'**\n');

%linear beam is B31
fprintf(fid,'*ELEMENT, TYPE=B21, ELSET=EL_fiber\n'); % For beam elements
fprintf(fid,'%8.0f,%8.0f,%8.0f\n',final_fibers');
fprintf(fid,'*ELEMENT, TYPE=CPS3, ELSET=EL_particle\n'); % For particle.  constant plane stress 3 node element
fprintf(fid,'%8.0f,%8.0f,%8.0f,%8.0f\n',particle_tri'); 

% --- SURFACES ---
fprintf(fid,'*SURFACE, NAME=P_Surf, TYPE=ELEMENT\n');
fprintf(fid,'EL_particle, \n');
fprintf(fid,'*SURFACE, NAME=N_Surf, TYPE=NODE\n');
fprintf(fid,'CONNECT, \n');

% --- MPC's ---
fprintf(fid,'**\n');
fprintf(fid,'** Multi Point Constraints\n');
fprintf(fid,'**\n');
fprintf(fid,'*TIE,NAME=INTERFACE, TYPE=NODE TO SURFACE, POSITION TOLERANCE = 0.1\n');
fprintf(fid,'N_Surf, P_Surf\n');

% % --- SECTION ASSIGNMENT AND MATERIAL PROPS ---

fprintf(fid,'**\n');
fprintf(fid,'** Section Assignment and Material Properties\n');
fprintf(fid,'**\n');

fprintf(fid,'** === === Beam Material Properties === ===\n');
fprintf(fid,'*BEAM GENERAL SECTION, ELSET=EL_fiber, DENSITY=%20.10E, SECTION=GENERAL\n',density); 
% No material given (it's defined here). Density is a required parameter.
% Poisson's ratio is optional; leaving blank sets it to zero.
fprintf(fid,'%20.10E,%20.10E,%20.10E,%20.10E,%20.10E\n',Area,MOI,0,MOI,J); % Area | I11 | I12 | I22 | J
fprintf(fid,'\n'); % blank line to accept default direction consines
fprintf(fid,'%20.10E,%20.10E\n',E,E/(2*(1+nu))); % Young's modulus | shear modulus | other entries unneeded
fprintf(fid,'*SECTION POINTS\n'); % Section points are needed to compute stress and strain within element
fprintf(fid,'%10.4f,%10.4f\n',0.5,0); % This gives x1 and x2 position of section points
fprintf(fid,'*DAMPING, ALPHA=%20.10E\n',0.007);


fprintf(fid,'** === === Particle Material Properties === ===\n');
fprintf(fid,'*MATERIAL, NAME=ThermalParticle\n');
fprintf(fid,'*ELASTIC, TYPE=ISOTROPIC\n');
fprintf(fid,'%20.10E,%20.10E,%20.10E\n',E*100,nu,0); %Modulus - Poisson - Initial Temp
fprintf(fid,'*EXPANSION, TYPE=ISO\n');
Alpha = 0.01; % 1 percent for 1 degree
fprintf(fid,'%20.10E,%20.10E\n',Alpha,0); %Exapnsion Coeffiction - Initial Temp
fprintf(fid,'*SOLID SECTION, ELSET=EL_particle, MATERIAL=ThermalParticle\n'); 

% 
% --- BOUNDARY CONDITIONS ---
fprintf(fid,'**\n');
fprintf(fid,'** Boundary Conditions\n');
fprintf(fid,'**\n');

% Fixed
fprintf(fid,'*BOUNDARY\n');
fprintf(fid,'NBNDRY, ENCASTRE\n');
fprintf(fid,'*BOUNDARY\n');
fprintf(fid,'INCLUSION, ENCASTRE\n');
fprintf(fid,'*INITIAL CONDITIONS, TYPE=TEMPERATURE\n');
fprintf(fid,'PARTICLE, 0\n');


%---Load Steps---
fprintf(fid,'**\n');
fprintf(fid,'** Load Steps\n');
fprintf(fid,'**\n');

% --- STEP 1 ---
fprintf(fid,'*STEP,NAME=STEP00%.3d, NLGEOM=YES, INC=300\n',1);
fprintf(fid,'*DYNAMIC, APPLICATION=QUASI-STATIC\n');
fprintf(fid,'%10.4f,%8.0f,%20.3e,%8.2f\n',[dt duration 0 0]); 

% Vary amplitude over time
%fprintf(fid,'*AMPLITUDE, NAME=LINEAR_AMPLITUDE_001, DEFINITION=TABULAR, TIME=STEP TIME, VALUE=RELATIVE\n');
fprintf(fid,'*AMPLITUDE, NAME=LINEAR_AMPLITUDE_001, DEFINITION=SMOOTH STEP, TIME=STEP TIME, VALUE=RELATIVE\n');
fprintf(fid,'%8.4f,%8.4f\n',[0 0]); % time 1 | load 1
fprintf(fid,'%8.4f,%8.4f\n',[ramp_time 1]); % time 2 | load 2


fprintf(fid,'*TEMPERATURE\n');
fprintf(fid,'PARTICLE, %8.0f\n',-Strain); 


% Outputs
fprintf(fid,'*OUTPUT, FIELD, NUMBER INTERVAL=5\n');
fprintf(fid,'*NODE OUTPUT\n');
fprintf(fid,'U, V, A, RF\n');
fprintf(fid,'*ELEMENT OUTPUT\n');
fprintf(fid,'LE, S, SK\n'); % SK gives section curvature change
fprintf(fid,'*OUTPUT, HISTORY, VARIABLE=ALL\n');

fprintf(fid,'*END STEP\n');


fclose('all');


end
