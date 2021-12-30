% This main script will extract the field output data from the ODB file in
% MATLAB and save the variables (output) in a MAT file
% It uses the extract_odb.m function, two python files
% written by Mainak Sarkar

% mention the odb file which is to be extracted:
odb_file_name = 'lattice_fiber_2020_07_07_sample_1_odb' ; % do not put extension, just mention name, if ODB is present in matlab current directory

% Specify the filename to save output as a MAT file
F = 'lattice_fiber_2020_07_07_sample_1_final_mat' ; % It is in the current directory of Matlab

%% extracting data from the ODB file
% Please go to the function extract_odb_file.m to specify what output variables (2 variables have been extracted, can be extended if 
% you need more variable extraction) you want and their components
% Please read the comments on the function file extract_odb_file to see how
% outputs are being saved in a text file.
read_table = 1 ; % 0, or 1; If you make this 0, then access output from three text files: 2 contain outputs and one contains frame time data
[U_comp_full, other_comp_full, frame_id_time] = extract_odb_file(odb_file_name, read_table) ;

%% saving output variables in a mat file
Path_mat = char(strcat(F,'.mat')) ;
save(Path_mat, 'U_comp_full', 'other_comp_full', 'frame_id_time') % output vars are in table form
