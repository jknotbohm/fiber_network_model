% Get Field Output From Specific Node Set
clear; clc; close all;

% Folder = uigetdir;
Folder='F:\Stephen\2019-Constant\BendingFull';

% Python Code and Arguments       
python_code = 'ODB_Elem.py';
ES = 'EL_fiber';
FO = 'S E' ;

copyfile(python_code, Folder);
copyfile('get_FO_request.py', Folder);
od=cd(Folder);
odb_list=dir('*.odb');
for i=1:length(odb_list)
    ODB=odb_list(i).name;

    CMD_Full = ['abaqus python ' python_code ' ' ODB ' ' ES ' ' FO];
    system(CMD_Full);

end
movefile(python_code, od);
movefile('get_FO_request.py', od);
cd(od);