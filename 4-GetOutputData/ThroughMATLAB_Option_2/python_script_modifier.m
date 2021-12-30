% extracts data from ODB files through Python scripts
% this function is written by Mainak Sarkar


function [] = python_script_modifier(F, test_no, field_out_var, index_field_out_var, fetch_time)


% Read text of the .py file into cell A
if fetch_time == 0
fid = fopen('C:\temp\read_odb_files.py','r');
else
fid = fopen('C:\temp\read_time_steps_frames.py','r');    
end
    
i = 1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end
fclose(fid);


% field_out_var = "U" ; % the field output variable in question
% index_field_out_var = "0" ; % count starts from zero (***NOTE***)
% test_no = "test2" ; % sequence it as per filed output variable considered

% F = 'sample' ;
% Change cell A
A{8} = strcat('odb_name=',"'",F,"'") ; % the name of the odb in the current directory of MATLAB
char(A{8}) ;
A{11} =  strcat('results_name=odb_name+',"'",test_no,"'") ;
char(A{11}) ;
if fetch_time == 0
A{36} = strcat('	odbSelectResults = x.fieldOutputs[',"'",field_out_var,"']"') ;
char(A{36}) ;
A{47} = strcat('		temp.append(s.data[',index_field_out_var,']) # the values in s.data are organized as follows: [S11,S22,S33,S12,S13,S23]') ;
char(A{47}) ;
end
% Write cell A into txt
if fetch_time == 0
fid = fopen('C:\temp\read_odb_files_modified.py', 'w');
else 
fid = fopen('C:\temp\read_time_steps_frames_modified.py', 'w');    
end
for i = 1:(numel(A)-1)
        fprintf(fid,'%s\n', A{i});
end