% Get Field Output From Specific Node Set
clear; clc; close all;

% Folder = uigetdir;
Folder_set = ["F:\Stephen\2019-Bending\2019-6-28 Linear High G3\K10-1";...
               "F:\Stephen\2019-Bending\2019-6-28 Linear High G3\K10-2";...
               "F:\Stephen\2019-Bending\2019-6-28 Linear High G3\K10-3";...
               "F:\Stephen\2019-Bending\2019-6-28 Linear High G3\K10-4";...
               "F:\Stephen\2019-Bending\2019-6-28 Linear High G3\K100";...
               "F:\Stephen\2019-Bending\2019-6-28 Linear High G3\K101";...
               "F:\Stephen\2019-Bending\2019-6-28 Linear High G3\K102";...
               ...
               ];
           
Folder_Main = 'F:\Stephen\2019-Bending\2019-2-12 Bending Full'
Folder_set = dir(Folder_Main);
Folder_set = Folder_set(3:end);
Folder_set= ["F:\Stephen\2019-Bending\2019-2-12 Bending Full";...
             "F:\Stephen\2019-Torsion\2019-2-1 Torsion Full"];
          
          
wbf = waitbar(0, 'Folders Completed')    ; 
wbj = waitbar(0, 'Jobs Completed');
for f = 1:length(Folder_set)
waitbar(f/length(Folder_set),wbf);
 Folder=char(string(Folder_set(f)));
% Folder = fullfile(Folder_set(f).folder,Folder_set(f).name);
% Folder = [Folder , '\Uniaxial Tension Tests'];
% Python Code and Arguments       
python_code = 'ODB_2_table.py';
NS = 'XPOS';
FO = 'U RF' ;

copyfile(python_code, Folder);
copyfile('get_FO_request.py', Folder);
od=cd(Folder);
odb_list=dir('*.odb');
for i=1:length(odb_list)
    waitbar(i/length(odb_list),wbj);
    ODB=odb_list(i).name;
%     Generate CMD & Run
%     if ~exist('results', 'var') 
%         CMD_Full = ['abaqus python get_FO_request.py ' ODB];
%         [~,results]=system(CMD_Full);
%         [FO,~] = listdlg('PromptString','Select Desired Outputs','ListString',strsplit(results,'\n'));
%         FO = num2str(FO);
%     end
    CMD_Full = ['abaqus python ' python_code ' ' ODB ' ' NS ' ' FO];
    system(CMD_Full);

end
movefile(python_code, od);
movefile('get_FO_request.py', od);
cd(od);
end