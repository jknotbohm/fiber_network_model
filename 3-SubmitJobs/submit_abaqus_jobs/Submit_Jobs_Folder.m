clear; close all; clc;
cd('..');
od=cd;
disp('Select Folder:')
new_Fold = uigetdir;
Folders{1} = new_Fold;
YN='y';
while true
    clc;
    disp('Current Folder List: ');
    for i = 1:length(Folders)
        disp(Folders{i});
    end
    if new_Fold == 0
        YN = input('Add Another Folder? (Y/N):   ','s');
    end
    if strcmpi(YN,'y')
        new_Fold = uigetdir;
        if new_Fold == 0 %GETDIR = Cancelled
            continue
        else
            Folders = [Folders; new_Fold]; %#ok<AGROW>
        end
    else
        break
    end
end

all_inp=[];
for i=1:length(Folders)
   cd(Folders{i});
   new_inp = dir('*.inp');
   all_inp = [all_inp;new_inp]; %#ok<AGROW>
end
cd(od);
h1 = figure; 
Initialise_INP_Pick(h1,all_inp)



% Functions 
function Initialise_INP_Pick(h,all_inp)
% h is figure handle
% Create INP List
W = h.Position(3); H = h.Position(4); 
LS.W=W; LS.H=H;
IL = uicontrol('Style','listbox','Position',[0.05*W 0.05*H 0.35*W 0.85*H],'String',{all_inp.name});
uicontrol('Style','text','Position',[0.05*W 0.9*H 0.35*W 0.05*H],...
           'String', ['Initial .inp files (',num2str(length(all_inp)),')']);
% Create Parameter Setup
AND_str = uicontrol('Style','edit','Position',[0.50*W 0.85*H 0.45*W 0.05*H]);
uicontrol('Style','text','Position',[0.50*W 0.90*H 0.45*W 0.05*H],...
           'String', 'Mandatory String Inclusion');
OR_str = uicontrol('Style','edit','Position',[0.50*W 0.70*H 0.45*W 0.05*H]);
uicontrol('Style','text','Position',[0.50*W 0.75*H 0.45*W 0.05*H],...
           'String', 'Must Include One');
% Create Update and Submit Buttons
uicontrol('Style','pushbutton','Position',[0.50*W 0.60*H 0.20*W 0.05*H],...
           'String', 'Update List','Callback',{@UpdateList,h});
uicontrol('Style','pushbutton','Position',[0.75*W 0.60*H 0.20*W 0.05*H],...
           'String', 'Submit Jobs','Callback',{@SubmitJobs,h});
% Save to Fig
LS.AND=AND_str;
LS.OR =OR_str;
LS.INP = all_inp;
guidata(h,LS)
UI_FinalList(LS.INP,H,W);
end


function UI_FinalList(inp_list,H,W)
uicontrol('Style','listbox','Position',[0.50*W 0.05*H 0.45*W 0.45*H],'String',{inp_list.name});
uicontrol('Style','text','Position',[0.50*W 0.50*H 0.45*W 0.05*H],...
           'String', ['Final .inp files (',num2str(length(inp_list)),')']);
end


function UpdateList(~,~,h)
    LS = guidata(h);
    AND = split(LS.AND.String);
    OR  = split(LS.OR.String);
    final_inp=LS.INP;
    for i=1:length(AND)
        AND_Del = ~contains({final_inp(:).name},AND(i));
        final_inp(AND_Del)=[];
    end
    OR_Del = ~contains({final_inp(:).name},OR);
    final_inp(OR_Del)=[];
    LS.INP=final_inp;
    guidata(h,LS);
    UI_FinalList(LS.INP,LS.H,LS.W);
end

function SubmitJobs(~,~,h)
    LS=guidata(h);
    close all;
    h2 = figure('Position',[580 558 760 420]);
    W = h2.Position(3); 
    H = h2.Position(4);
    %Initialize Table
    InputName = {LS.INP.name}';
    [Start,End,RunTime] = deal(zeros(length(InputName),1));
    Complete = string(false(length(InputName),1)); Complete(:)='Pending';
    T = table(InputName,Start,End,RunTime,Complete);
    uit=uitable(h2,'Data',cellstr([InputName,Start,End,RunTime,Complete]),'Position',[0.05*W 0.05*H 0.90*W 0.85*H]);
    uit.ColumnName = {'Input Name'; 'Start Time'; 'End Time'; 'Total'; 'Status'};
    datetime.setDefaultFormats('default','hh:mm a');
    for i=1:length(InputName)
%         figure(h2)
%         waitbar(i/length(InputName),InputName(i));
        tic;
        uit.Data(i,2)=cellstr(datetime);
        uit.Data(i,5)={'Running'};
        pause(1)
        % Run Job
        folder=LS.INP(i).folder;
        od=cd(folder);
        title(folder)
        inp_name = LS.INP(i).name(1:end-4);
        if ~isfile([inp_name,'.odb'])
        command=['cd ',folder];
        system(command);
        command = ['abaqus job=' ,inp_name, ' interactive ask_delete=off'];
        system(command);
        end
        % Fill Remainder
        uit.Data(i,3)=cellstr(datetime);
        t=toc;
        uit.Data(i,4) = {[num2str(round(t/60)),'m']};
        uit.Data(i,5) = {'Completed'};
        cd(od)
    end
end