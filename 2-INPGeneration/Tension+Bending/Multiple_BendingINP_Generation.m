close all; clear; clc;   

Date = char(datetime('now','Format','yyyy-MM-dd'));
SimTitle = 'Bending Tests';
Folder = [Date ' ' SimTitle];
clear Date SimTitle
mkdir(Folder)

theta = .1*pi/180; %Given in radians
strain = .01; %Uniaxial Strain

L_set = 30;
R_set = 5:10;
index_set = 1001:1003;
    
%Create Cylinders and INP
for i = (index_set)
    file_loc=Call_CYL_INP_Generation(Folder, L,R,i);
    save_loc = [file_loc '\Uniaxial Tension Tests'];
    Static_2D_PureBending_INP(file_loc, theta);
    Static_Tension_INP_Gen(file_loc,save_loc,strain)
end
