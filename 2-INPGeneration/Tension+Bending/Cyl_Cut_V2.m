function [ final_nodes,final_fibers] = Cyl_Cut_V2( final_nodes,final_fibers,L,R)
%Rewritten Version Of Cyl Cut to improve efficiency

final_rad = sqrt(final_nodes(:,3).^2 + final_nodes(:,4).^2);
% Indexed Arrays
nodes_outside_L = abs(final_nodes(:,2)) > L/2;
nodes_outside_R = abs(final_rad) > R;
nodes_outside_OR=  or(nodes_outside_L,nodes_outside_R);
% nodes_inside_OR= or(~nodes_outside_L,~nodes_outside_R); %Useless Array


% Move all nodes outside the boundary to the boundary.
% This makes a more uniform domain witout interfering with deletion
% Move outside L Nodes
final_nodes(nodes_outside_L,2) = ...
   final_nodes(nodes_outside_L,2)./abs(final_nodes(nodes_outside_L,2))*L/2;
% Move outside R Nodes
final_nodes(nodes_outside_R,[3 4]) = ...
   final_nodes(nodes_outside_R,[3 4])./final_rad(nodes_outside_R)*R;


% Determine if Fibers are Completely out, completely in, or border
Nodes_A = final_fibers(:,2);
Nodes_B = final_fibers(:,3);
Comp_Out = nodes_outside_OR(Nodes_A) & nodes_outside_OR(Nodes_B);
% If not Comp_Out, then it will not be deleted.  Therefore these are not
% needed
% Comp_In = nodes_inside_OR_ind(Nodes_A) & nodes_inside_OR_ind(Nodes_B);
% Border = abs(nodes_outside_OR_ind(Nodes_A) - nodes_outside_OR_ind(Nodes_B));
%     % If both in, 0-0=0.   If both out 1-1=0.  Else =1

%RemoveDanglingFib
n_Unique = 1;
while n_Unique ~= 0
         edges = 0.5:1:(max(final_nodes(:,1))+0.5);
         N = histcounts(final_fibers(:,2:3),edges);
         single_nodes = find(N==1);
         f2d1 = ismember(final_fibers(:,2),single_nodes);
         f2d2 = ismember(final_fibers(:,3),single_nodes);
         f2d = or(f2d1,f2d2);
         n_Unique = sum(f2d);
         if n_Unique == 0
             continue
         end
         final_fibers(f2d,:)=[NaN, NaN, NaN]; %#ok<SAGROW>
end


%Remove Unused Fibers
final_fibers(Comp_Out,:)=NaN;
used_nodes = unique([final_fibers(:,2),final_fibers(:,3)]);
unused_nodes = setdiff(final_nodes(:,1),used_nodes);
final_nodes(unused_nodes,:)=NaN;

end