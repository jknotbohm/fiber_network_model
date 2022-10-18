
%% This function determines the number of nodesalong x and y directions.
% INPUT: density and domain dimensions of the network
% OUTPUT: numbers and spacings of nodes
% Written by Mainak Sarkar, University of Wisonsin-Madison

function [Nx, Ny, dNx, dNy] = find_node_numbers_2Dlattice(rho_nodal, lx, ly)
syms Nx Ny real
dNx = lx/(Nx-1);            % Spacing in x direction
dNy = dNx*sqrt(3)/2;        % Spacing in y direction (not equal for a hex connectivity)
Ny = (ly/dNy)+1;            % Number of nodes in y direction

fun = ( (Nx * Ny) / (lx * ly) ) - rho_nodal ;

sol = solve([fun == 0], [Nx]) ;

if sol(1) > 0
    Nx = round(double(sol(1))) ;
else 
    Nx = round(double(sol(2))) ;
end

dNx = lx/(Nx-1);            % Spacing in x direction
dNy = dNx*sqrt(3)/2;        % Spacing in y direction (not equal for a hex connectivity)

Ny = round(ly/dNy)+1;             % Number of nodes in y direction