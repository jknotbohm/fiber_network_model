% This script generates the text file containing field output variable
% values in a comma separated txt file, can be visualized on Notepad ++ and
% excel, for all frames (rows) and for all the nodes (columns)
% written by Mainak Sarkar


function [U_comp_full, other_comp_full, frame_id_time] = extract_odb_file(F, read_table)

% U data for all frames and for all nodes:
test_no = "test11" ; % sequence it as per filed output variable considered
field_out_var = "U" ; % the field output variable in question
index_field_out_var = "0" ; % count starts from zero (***NOTE***). Eg. U1 will be "0", U2 will be "1".
step_no = "0" ; % step 1 
python_script_modifier(F, test_no, step_no, field_out_var, index_field_out_var, 0) 
system('abaqus python C:\temp\read_odb_files_modified.py');
Path1 = char(strcat(F,test_no,'.txt')) ;
if read_table == 1
 U_comp_full = load(Path1) ;
else
 U_comp_full = NaN ;
end

% Cauchy stress data for all frames and for all nodes (integration points):
test_no = "test12" ; % sequence it as per filed output variable considered
field_out_var = "RF" ; % the field output variable in question, RF
index_field_out_var = "0" ; % count starts from zero (***NOTE***). Eg. mag. RF1 will be "0", RF2 will be "1".
step_no = "0" ; % step 1 
python_script_modifier(F, test_no, step_no, field_out_var, index_field_out_var, 0)
system('abaqus python C:\temp\read_odb_files_modified.py');
Path2 = char(strcat(F,test_no,'.txt')) ;
if read_table == 1
 other_comp_full = load(Path2) ;
 else
 other_comp_full = NaN ;
end

% Following lines generate a text file 'frame_time'
step_no = "0" ; % step 1 
python_script_modifier(F, test_no, step_no, field_out_var, index_field_out_var, 1) 


Path = [F,'_frame_times.txt'] ;
diary(Path)
system('abaqus python C:\temp\read_time_steps_frames_modified.py');
diary off

if read_table == 1
 frame_id_time = load(Path);
else
 frame_id_time =  NaN ;
end

return


