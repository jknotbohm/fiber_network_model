%% Fiber Network Creation Main Script;
%
%   Functions used by this script:
%       -Initial_Network_Generation
%       -Network_Optimization
%       -Branch_Optimization
%
% This script generates a fiber network that gets saved as a .mat file
%
%
% Written by Peter Grimmer
% Modified by Stephen Tyznik to increase efficiency 
% Modified to form 2D network directly by Mainak Sarkar

clear; close all; clc;


for File_Index_Inp = 1001:1003
    %Set RNG Seed
    rng(File_Index_Inp)      ;
     
    % Basic parameters:
    l_fiber = 0.1         ;    % desired fiber length, not the initial avg fiber length
    lx = 25                ;
    ly = 25               ;    % Domain dimensions
     rho_base = 6   ;    % This is the standard density used for most tests, default 0.6, 0.021, 1 is ok; >1 causes trouble in INP generators
     rho_scale = 1.6      ;    % node density: nodes/unit volume, default 1, 0.005
     rho = rho_base*rho_scale ;
    % rho = 0.6 ;
    
    % Filename
    fn = 'Network_';
    File_Index = num2str(File_Index_Inp);
    ext = '.mat';
    Network_Variables = ['2DNetwork_',num2str(lx),'x' num2str(ly) '_Seed_' File_Index];

    % Also going to save the network before branching optimization:
    ub = '_unbranched';
    Network_variables_unbranched = [fn File_Index ub ext];

    % Set up paths to save files to:
    curdir = pwd;
    mkdir(Network_Variables)% make a directory for the file index
    path = [curdir '\' Network_Variables '\' ];


    %% Initial Network Generation

    [nodes, fibers, fiberenergy, total_energy, N1, N2, N_boundary_nodes,fiberlengths,valencycheck] = Initial_Network_Generation_2D(rho,lx,ly,l_fiber,200);

    % nodes: N x 3, nodal coordinates
    % fibers: N_fibers x 2, fiber node numbers
    % total_energy: scalar measure of fiber length. Network_Optimization will
    %      minimize
    N = length(nodes);
    N_fibers = length(fibers);
    str = ['There are ' num2str(N) ' nodes and ' num2str(N_fibers) ' fibers'];
    disp(str);
    str = ['The total Initial Fiber Length Energy is ' num2str(total_energy)];
    disp(str);


    %% Fiber Length Optimization
    N_anneal=15;
    N_optimize = 5;
    stepsize = linspace(1, .1, N_optimize);
    fraction_to_try_swap = linspace(.3,.02,N_optimize);

    for j=1:N_optimize
        swap_skip_energy = 3*(median(fiberlengths)-l_fiber)^2;
        str=['Swap skip energy threshold = ' num2str(swap_skip_energy)];
        disp(str);
        [nodes, fibers, fiberenergy, fiberlengths, total_energy, N1,N2,N_interior] = Network_Optimization_2D(fraction_to_try_swap(j),N,nodes,fibers,N_anneal,lx,ly,l_fiber,fiberlengths,fiberenergy,N1,N2, N_boundary_nodes,stepsize(j),swap_skip_energy);
        percent_done = 100*j/N_optimize;
        str=['percent optimized = ' num2str(percent_done)];
        disp(str);
        str='_________________________________________';
        disp(str);
        save([path Network_variables_unbranched]);  % Save each loop iteration overwriting the previous file each time
    end

    figure
    histogram(fiberlengths);
    title('Fiber Length Distribution after Optimization')


    % Branching Optimization
    N_branching_optimize = 4;
    [nodes, fibers, nodal_branching_energy, total_branching_energy_init, total_branching_energy_final] = Branch_Optimization_2D(N_branching_optimize, nodes, fibers, N);


    save([path Network_Variables]);

end
