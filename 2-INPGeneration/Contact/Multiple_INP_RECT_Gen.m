close all; clear; clc;   

Date = char(datetime('now','Format','yyyy-MM-dd'));
SimTitle = 'Rigid Inclusion Displacement Tests';
Dir = [Date ' ' SimTitle];
clear Date SimTitle
mkdir(Dir)

R=50:50:450;
index_set = 1001:1004;
%Create Cylinders and INP
    for i = (index_set)
        for R_IND=R
            LY=2*R_IND;
            LX=4*R_IND;
            w = R_IND/50;
            file_loc=Call_SEMICYL_INP_Generation(Folder, LX, LY, i);
            Static_2D_Cyl_Ind_INP(file_loc, w, R_IND)
        end 
    end
