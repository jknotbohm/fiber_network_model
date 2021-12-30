close all; clear; clc; 


% load(['C:\Users\msarkar3\Desktop\MATLAB Path\3D_printing_2D\3DNetwork_25x25x2_Seed_1003\3DNetwork_25x25x2_Seed_1003.mat'])



Date = char(datetime('now','Format','yyyy-MM-dd'));
SimTitle = 'utx_v4_23Feb';
Folder = [Date ' ' SimTitle];
clear Date SimTitle
mkdir(Folder)

theta = .1*pi/180 ; %Given in radians
strain = 0.3 ; %Uniaxial Strain / simple shear strain (depends on case)

L_set = 350 ;
R_set = 430 ;
index_set = 1003;
 [nodes, fibers, lx, ly, lz] = auto_crop_tool(0.3) ; 
% dangling fibers are removed. 
curr_density = size(nodes,1)/(lx*ly*lz) ;
disp(['post-crop density = ', num2str(curr_density)])

% Create Cylinders and INP
for i = (index_set)
     [file_loc, l_fiber] = Call_CYL_INP_Generation(Folder, L_set, R_set, i);
     save_loc = [file_loc '\utx_v4_23Feb'];
%     Static_2D_PureBending_INP(file_loc, theta);
%    [final_nodes, final_fibers] = Static_Tension_INP_Gen(file_loc,strain, l_fiber) ;
% [final_nodes, final_fibers] = Static_Tension_INP_Gen_STEP(file_loc,strain, l_fiber)  ;
[final_nodes, final_fibers] = Static_Tension_INP_Gen_STEP_FxyCz(file_loc,strain, l_fiber)  ;
% [final_nodes, final_fibers] = Static_Tension_INP_Gen_STEP_X(file_loc,strain, l_fiber)  ;
%      [final_nodes, final_fibers] = Simple_Shear_INP_Gen(file_loc,strain, l_fiber) ; % for simple shear
% [final_nodes, final_fibers] = Simple_Shear_INP_Gen_STEP(file_loc,strain, l_fiber) ;
% [final_nodes, final_fibers] = triaxial_INP_Gen(file_loc,strain, l_fiber) ; % triaxial
% [final_nodes, final_fibers] = Triaxial_INP_Gen_STEP(file_loc,strain, l_fiber) ;
% [final_nodes, final_fibers] = Simple_Shear_INP_Gen_STEP_nbul(file_loc,strain, l_fiber) ;
end


index_inp = 1003 ;
save(['C:\Users\msarkar3\Desktop\MATLAB_Path\3D_printing_2D\',...
    '3D_d0p01f30_Network_350X350x5_Seed_',num2str(index_inp),'\'...
    '3D_d0p01f30_Network_350x350x5_Seed_cropped_',num2str(index_inp),'.mat'],'final_nodes','final_fibers', '-append') % #ok<LOAD>
