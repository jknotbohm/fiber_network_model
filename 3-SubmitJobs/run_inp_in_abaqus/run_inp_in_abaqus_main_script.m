% This main script runs any INP file in Abaqus, stores the result in the
% MATLAB current directory as an ODB file 
% written by Mainak Sarkar

% load the path of INP file to be run on Abaqus
file = 'C:\Users\msarkar3\Desktop\MATLAB_Path\3D_printing_2D\lattice_fiber_2020_07_07_sample_1_inp\strain_0.06.inp' ; 

% name the ODB file name
odb_file_name = 'lattice_fiber_2020_07_07_sample_1_odb' ; % This ODB file will be stored in the current folder (if only filename is mentioned, not the full path) of MATLAB after creation in the Abaqus CAE post-simulation


%% MATLAB-Abaqus CAE interfacing: Run the INP file in Abaqus and generate odb file 
% in the MATLAB current folder
run_inp_abaqus(file, odb_file_name)

