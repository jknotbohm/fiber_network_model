function [final_nodes,final_fibers] = Rect_Cut_V2(final_nodes,final_fibers,lx,ly,lz)

%Crop_to_Sphere: Crops the box domain into a spherical domain
%   Search for all nodes and fibers outside of a given outer and inner
%   radius

nodes_outside_X = abs(final_nodes(:,2)) > lx/2;
nodes_outside_Y = abs(final_nodes(:,3)) > ly/2;
nodes_outside_Z = abs(final_nodes(:,4)) > lz/2;
nodes_outside_OR = or(nodes_outside_Z,or(nodes_outside_X,nodes_outside_Y));

% Move all nodes outside the boundary to the boundary.
% This makes a more uniform domain witout interfering with deletion
% Move outside Lx Nodes
final_nodes(nodes_outside_X,2) = ...
   final_nodes(nodes_outside_X,2)./abs(final_nodes(nodes_outside_X,2))*lx/2;
% Move outside Ly Nodes
final_nodes(nodes_outside_Y,3) = ...
   final_nodes(nodes_outside_Y,3)./abs(final_nodes(nodes_outside_Y,3))*ly/2;
% Move outside Lz Nodes
final_nodes(nodes_outside_Z,4) = ...
   final_nodes(nodes_outside_Z,4)./abs(final_nodes(nodes_outside_Z,4))*lz/2;


% Determine if Fibers are Completely out, completely in, or border
Nodes_A = final_fibers(:,2);
Nodes_B = final_fibers(:,3);
Comp_Out = nodes_outside_OR(Nodes_A) & nodes_outside_OR(Nodes_B);
% If not Comp_Out, then it will not be deleted.  Therefore these are not
% needed
% Comp_In = nodes_inside_OR_ind(Nodes_A) & nodes_inside_OR_ind(Nodes_B);
% Border = abs(nodes_outside_OR_ind(Nodes_A) - nodes_outside_OR_ind(Nodes_B));
%     % If both in, 0-0=0.   If both out 1-1=0.  Else =1

%Remove Unused Fibers
final_fibers(Comp_Out,:)=NaN;
used_nodes = unique([final_fibers(:,2),final_fibers(:,3)]);
unused_nodes = setdiff(final_nodes(:,1),used_nodes);
final_nodes(unused_nodes,:)=NaN;
end
