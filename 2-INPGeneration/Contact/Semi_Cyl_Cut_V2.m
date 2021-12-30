function [final_nodes,final_fibers,R_network] = Semi_Cyl_Cut_V2(final_nodes,final_fibers,LX,LY)
% Takes rectangular domain and cuts upper edge to a circular domain
% Cut will be made such that Y(x=0) = LY/2, Y(x=+-LX/2)=-LY/2

% %% DEBUG CODE
% load('F:\Stephen\2019-Cyl Punch\2019-2-7 CONTACT Hertzian\400_x_200_Plate_1003.mat','final_nodes','final_fibers');
% LX=100;
% LY=50;

%% Solve for Arc to Cut
% Solve for circular Radius, and Center
Bottom = -LY/2 +0.75*LY;
syms yo R theta
eq1 =  LY/2 == R*sin(pi/2) + yo; %#ok<*NODEF>
eq2 =  LX/2 == R*cos(theta);
eq3 = Bottom == R*sin(theta) + yo;
[yo, R_network, theta] = solve([eq1,eq2,eq3],[yo,R,theta]);
yo=double(yo); R_network=double(R_network); theta=double(theta);
tt=linspace(theta,pi-theta);
xx=R_network*cos(tt); yy=R_network*sin(tt)+yo;


%% First: Out of Bounds Nodes
nodal_R = sqrt(final_nodes(:,2).^2 + (final_nodes(:,3)-yo).^2);
nodes_outside_X = abs(final_nodes(:,2)) > LX/2;
nodes_outside_Y = final_nodes(:,3) < -LY/2;
Nodes_Outside_CYL = nodal_R > R_network;
Nodes_Outside_FULL = or(Nodes_Outside_CYL,or(nodes_outside_X,nodes_outside_Y));
% Move all nodes outside the boundary to the boundary.
% This makes a more uniform domain witout interfering with deletion
% Move outside Lx Nodes

%% Cut Network To Be within New Domain
old_nodes = final_nodes;
final_nodes(nodes_outside_X,2) = ...
   final_nodes(nodes_outside_X,2)./abs(final_nodes(nodes_outside_X,2))*LX/2;
% Move outside Ly Nodes
final_nodes(nodes_outside_Y,3) = ...
   final_nodes(nodes_outside_Y,3)./abs(final_nodes(nodes_outside_Y,3))*LY/2;
% Move Outside R Nodes - X
final_nodes(Nodes_Outside_CYL,2) = ...
   final_nodes(Nodes_Outside_CYL,2)./nodal_R(Nodes_Outside_CYL)*R_network;
% Move Outside R Nodes - Y
final_nodes(Nodes_Outside_CYL,3) = ...
   (final_nodes(Nodes_Outside_CYL,3)-yo)./nodal_R(Nodes_Outside_CYL)*R_network+yo;


% Determine if Fibers are Completely out, completely in, or border
Nodes_A = final_fibers(:,2); 
Nodes_B = final_fibers(:,3); 
Comp_Out = Nodes_Outside_FULL(Nodes_A) & Nodes_Outside_FULL(Nodes_B);
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

%% Plot for sanity
% close all
% plot(old_nodes(:,2),old_nodes(:,3),'bo',final_nodes(:,2),final_nodes(:,3),'ro')
% figure
% plot(old_nodes(Nodes_Outside_CYL,2),old_nodes(Nodes_Outside_CYL,3),...
%     'bo',final_nodes(Nodes_Outside_CYL,2),final_nodes(Nodes_Outside_CYL,3),'ro'),...
%     xx,yy,'k--',0,yo,'k+'); axis equal

end