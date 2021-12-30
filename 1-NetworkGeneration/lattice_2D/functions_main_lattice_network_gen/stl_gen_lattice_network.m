% this function file generates the stl file from the given network
% written by Mainak Sarkar

function [] = stl_gen_lattice_network(el_set_final, nodes_set_final, d)


N = 100 ;  % cylindrical beam as fibers, do not change it.
R = d/2 ; % c/s radius of fibers

% make the folder that will contain the STL files of fibers and grips
% mkdir 100by100_cropped_NN_4
P = 'C:\Users\msarkar3\Desktop\MATLAB_Path\3D_printing_2D';
F = 'network_lattice_2020_07_20_trial3';
mkdir(P,F)
new = fullfile(P,F) ;


% to run this functtion file prior to adding midnodes to each fibers in
% main lattice generator script.
t = nodes_set_final(1:end,:) ;
[N_nodes, columns] = size(t) ;
t2 = [el_set_final(:,1), el_set_final(:,2), el_set_final(:,3)] ;
N_number = t(:,1) ;
x = t(:,2) ;
y = t(:,3) ; 

% generating the stl file
nn = 0 ;
figure(1)
counter = 0 ;
counter1 = 0 ;
for i = 1:size(t2,1)
   for j = 1:N_nodes
   if t2(i,2)==t(j,1)
       counter = counter + 1 ;
       X(counter,:) = [t(j,2) t(j,3)] ;
   end
   end
   for j = 1:N_nodes
   if t2(i,3)==t(j,1) 
       counter1 = counter1 + 1 ;
       Y(counter1,:) = [t(j,2) t(j,3)] ;
   end
   end
   nn = nn + 1 ; % calculates the number of fibers
[Xa,Ya,Za] = cylinder2P(R, N, [X(counter,:) 0], [Y(counter1,:) 0]) ;
% C{i} = [Xa, Ya, Za] ;
filename_reference = [new,'\fiber_no_',num2str(nn),'.stl'] ;
surf2stl(filename_reference,Xa,Ya,Za)
hold on 
end 

view(2)
daspect([1 1 1])
pbaspect([1 1 1])
axis off
 

%% generate the python script to run in blender python console to obtain final stl file of the network:
folder = ['C:/Users/msarkar3/Desktop/MATLAB_Path/3D_printing_2D/',F] ;
fiber_nos = nn 
final_stl_fileneme = 'merged_stl_file' ; % only finename / no extension here

blender_python_script_gen(fiber_nos, folder, final_stl_fileneme) ;

 