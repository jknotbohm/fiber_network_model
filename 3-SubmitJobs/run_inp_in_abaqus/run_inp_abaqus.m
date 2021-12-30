% running INP file from MATLAB in Abaqus CAE, and subsequently generating
% an .odb file which is output database. The .odb file is saved in the
% current folder. 
% written by Mainak Sarkar


function [] = run_inp_abaqus(inp_file_name, odb_file_name)
% two commented lines below are to be defined in the main file:
% inp_file_name = 'C:\Users\msarkar3\Desktop\MATLAB_Path\3D_printing_2D\lattice_fiber\strain_0.005.inp' ;
% odb_file_name = 'filename5' ;
s = ['abaqus job=',odb_file_name,' ','inp=',inp_file_name ,' ', 'interactive'] ;
dos(s)  % running the job in Abaqus CAE 

%-------------------------------------------------------------------------------------------------------------------------------------
% This is the source / format:
% dos('abaqus job=filename2 inp=C:\Users\msarkar3\Desktop\MATLAB_Path\3D_printing_2D\lattice_fiber\strain_0.005.inp interactive')


