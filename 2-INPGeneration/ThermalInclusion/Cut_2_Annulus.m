function [File_Path]=Cut_2_Annulus(Dest_Fold,Net_Path,Inner_Rad,Outer_Rad)

% Load in Netork File
load(Net_Path);
if Outer_Rad == Inf
    Outer_Rad=min(lx,ly)/2; %If undefined, make largest possible circular domain
end
Seed = File_Index_Inp;
% Define filename and save location
File_Name = ['Annulus_' num2str(Outer_Rad) 'x' num2str(Inner_Rad) '_Index_' Net_Path(end-7:end)];
File_Path = fullfile(Dest_Fold,File_Name);

% Remove any duplicated fibers
fibers=sort(fibers,2); %Put smallest index first for comparison
setOrder='stable';
[~,unique_fibers,~] = unique(fibers,'rows',setOrder);
fibers = fibers(unique_fibers,:);
clear setOrder unique_fibers

% Create 'Final' Arrays
final_nodes = [(1:length(nodes))' nodes];
final_fibers = [(1:length(fibers))' fibers];

% Cut Outer Border
[final_nodes,final_fibers] = cut_outer_bound(final_nodes,final_fibers,Outer_Rad);

%Cut Inner Border
[final_nodes,final_fibers] = cut_inner_bound(final_nodes,final_fibers,Inner_Rad);

%Save Final Network
save(File_Path,'final_nodes','final_fibers','Seed')

end

function [final_nodes,final_fibers] = cut_outer_bound(final_nodes,final_fibers,Outer_Rad)
%Find Nodes Within Boundsary
Nodes_Rad = sqrt(final_nodes(:,2).^2 + final_nodes(:,3).^2);
Nodes_Out = Nodes_Rad > Outer_Rad;
Nodes_Out_Label = final_nodes(Nodes_Out,1);

%Move Nodes To Boundary
final_nodes(Nodes_Out,2:3) = final_nodes(Nodes_Out,2:3) ./ Nodes_Rad(Nodes_Out) * Outer_Rad;

%Remove Completely Out Fibers
Comp_Out = ismember(final_fibers(:,2),Nodes_Out_Label) & ismember(final_fibers(:,3),Nodes_Out_Label);
final_fibers(Comp_Out,:) = [];

%Rename Arrays
[final_nodes,final_fibers] = rename_arrays(final_nodes,final_fibers);
end

function [final_nodes,final_fibers] = cut_inner_bound(final_nodes,final_fibers,Inner_Rad)
%Find Nodes Within Boundsary
Nodes_Rad = sqrt(final_nodes(:,2).^2 + final_nodes(:,3).^2);
Nodes_Out = Nodes_Rad < Inner_Rad;
Nodes_Out_Label = final_nodes(Nodes_Out,1);

%Find Intersect Fibers
Int = xor(ismember(final_fibers(:,2),Nodes_Out_Label),ismember(final_fibers(:,3),Nodes_Out_Label));
Int_Node_log = [0*Int, Int, Int] .* ismember(final_fibers,Nodes_Out_Label);
Int_Node_Ind = find(Int_Node_log); %Used for renaming final fibers array later;
[Int_row,Int_col_1] = ind2sub(size(final_fibers),Int_Node_Ind);
    Int_col_2 = -1*(Int_col_1-3)+2; %Switch all 2->3 and 3->2
%Find the out of bounds (1) and in bounds (2) nodes of the
N1 = final_fibers(sub2ind(size(final_fibers),Int_row,Int_col_1)); 
X1 = final_nodes(N1,2); Y1 = final_nodes(N1,3);
N2 = final_fibers(sub2ind(size(final_fibers),Int_row,Int_col_2)); 
X2 = final_nodes(N2,2); Y2 = final_nodes(N2,3);
%Find Line Equations
dy = (Y2-Y1);
dx = (X2-X1);
%solve with Quadratic formula
a = dy.^2 + dx.^2;
b = 2*Y1.*dy + 2*X1.*dx;
c = Y1.^2 + X1.^2 - Inner_Rad^2;
sol1 = (-b + sqrt(b.^2-4*a.*c))./(2*a);
sol2 = (-b - sqrt(b.^2-4*a.*c))./(2*a);
scale = max([sol1,sol2],[],2); %positive indicates going from node 1 to node 2.  should always be max
%New Positions and Labels
Xn = X1+dx.*scale;
Yn = Y1+dy.*scale;
Zn = zeros(size(Xn));
Ln = final_nodes(end,1) + (1:length(Xn))'; 
% add to final nodes
final_nodes = [final_nodes; [Ln,Xn,Yn,Zn] ];
% reaname fiber array
final_fibers(Int_Node_Ind) = Ln;

%Remove Completely Out Fibers
Comp_Out = ismember(final_fibers(:,2),Nodes_Out_Label) & ismember(final_fibers(:,3),Nodes_Out_Label);
final_fibers(Comp_Out,:) = [];

%Rename Arrays
[final_nodes,final_fibers] = rename_arrays(final_nodes,final_fibers);
end

function [final_nodes,final_fibers] = rename_arrays(final_nodes,final_fibers)
%Rename Fiber Labels
final_fibers(:,1) = (1:length(final_fibers))';

%Find Unique Nodes
Used_Nodes =  unique(final_fibers(:,2:3));
New_Nodes = (1:length(Used_Nodes))';

%Remove unused nodes and relabel
final_nodes(~ismember(final_nodes(:,1),Used_Nodes),:)=[];
final_nodes(:,1) = New_Nodes;

%Relabel Fiber Indexes
[~,Node_Ind]=ismember(final_fibers(:,2:3),Used_Nodes);
final_fibers(:,2:3) = New_Nodes(Node_Ind);
end