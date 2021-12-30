% this function generates a python script which has to be run on the
% Blender-Python console to obtain the final STL file for the entire
% network.
% PYTHON SCRIPT IS SAVED IN THE CURRENT FOLDER OF MATLAB.
% written by Mainak

function [] = blender_python_script_gen(fiber_nos, folder, final_stl_fileneme)
%{
folder = ['C:/Users/msarkar3/Desktop/MATLAB_Path/3D_printing_2D/',F] ;
fiber_nos = 300 ;
final_stl_fileneme = 'merged_stl_file' ; % only finename / no extension here
%}

clear A

A{1} = ['import os'] ;
A{2} = ['import bpy'] ;

A{3} = ['bpy.ops.object.select_all(action=','"','SELECT")'] ;
A{4} = ['bpy.ops.object.delete()'] ;

for gh = 1 : fiber_nos
str = join(strcat(folder,'/fiber_no_', char(string(gh)), '.stl')) ;
f = char(str) ;
A{5+gh} = ['bpy.ops.import_mesh.stl(filepath="',f ,'",','filter_glob="*.stl",files=[{"name": "fiber_no_',char(string(gh)),'.stl"}],directory="',folder,'/")'] ;
char(A{5+gh}) ;
end

A{5 + gh + 1} = join(["bpy.ops.object.select_all(action=","'SELECT')"]) ;
char(A{5 + gh + 1}) ;

A{5 + gh + 2} = join("bpy.ops.object.join()") ;
char(A{5 + gh + 2}) ;

path_stl = [folder, '/', final_stl_fileneme] ;
addem_str = char(["'OFF', axis_forward='Y', axis_up='Z')"]) ;
addem_str= addem_str(~isspace(addem_str));
A{5 + gh + 3} = ['bpy.ops.export_mesh.stl(filepath="',path_stl,'.stl", check_existing=True, filter_glob="', final_stl_fileneme,'.stl", use_selection=False, global_scale=1.0, use_scene_unit=False, ascii=False, use_mesh_modifiers=True, batch_mode=',char(addem_str)]  ; 


% Write cell A into txt
fclose('all') ;
if exist('C:\Users\msarkar3\Desktop\MATLAB_Path\3D_printing_2D\model_gen_script_blender_modified.py', 'file')==2
  delete('C:\Users\msarkar3\Desktop\MATLAB_Path\3D_printing_2D\model_gen_script_blender_modified.py');
end
fid = fopen('C:\Users\msarkar3\Desktop\MATLAB_Path\3D_printing_2D\model_gen_script_blender_modified.py', 'w');

for i = 1:(numel(A))
    fprintf(fid,'%s\n', A{i});
end

