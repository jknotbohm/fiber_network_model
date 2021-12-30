function [file_loc] = Call_RECT_INP_Generation(Folder, LX,LY,index_inp)
%% Network Modifications for INP Generation: Cylindrical domain for tension or torsion

% Crops a 3D cubic network into a cylidner and creates sets of nodes for
% BCs


load(['C:\Users\Stephen\Desktop\Abaqus\AbaqusMatlabCodes\Network Generation\',...
    'Tyz_1000x500x1_Plate_Rho_1x_',num2str(index_inp),'.mat\',...
    'Tyz_1000x500x1_Plate_Rho_1x_',num2str(index_inp),'.mat'])



%N_nodes_final should be N_nodes_initial*(pi*r^2)*height;

% Precondition the network for making the INP file. Run this then run the
% INP generation script directly after without clearing, or load the
% cropped network mat file before running INP generation script

% The network mat file has lx, ly, lz which are dimensions of the brick
% domain. Usually I make cubic domains so lx=ly=lz. This will crop a
% domain into a cylinder with z being the axial direction. THe radius of
% the cylinder should be prescribed in terms of lx or ly. If lx=ly, the
% biggest cylinder you would make would be radius=lx/2 or ly/2, and height
% lz

height=lz;

%% Delete any mirrored duplicate fibers that made it through the prior processes:

% Due to floating point node location values we can identify duplicate fibers by
% computing the absolute values of the fiber orientation unit vectors.

cropped_network = [num2str(LX) ,'_x_', num2str(LY),'_Plate_' ,num2str(index_inp)];

path=[cd,'\',Folder];
Orig_Dir=cd(path);
%  if exist([cropped_network,'.mat'], 'file') == 2
%      load([cropped_network,'.mat']);
%      path=[cd,'\'];
%      file_loc=[path,cropped_network];
%      save(cropped_network)
%      cd(Orig_Dir)
%      return
%  end
cd(Orig_Dir);

UV = zeros(length(fibers),3);

for j=1:length(fibers)
    
    uv = nodes(fibers(j,1),:) - nodes(fibers(j,2),:);
    L = sqrt(sum(uv.^2));
    uv = uv/L;
    UV(j,:) = abs(uv);  % Bc fibers could be oriented opposite but still the same

end

setOrder='stable';
[~,unique_fibers,~] = unique(UV,'rows',setOrder);

N_duplicate_fibers = length(fibers)-length(unique_fibers);

fibers = fibers(unique_fibers,:);

str = ['The number of duplicate fibers deleted was ' num2str(N_duplicate_fibers)];
disp(str);

N_fibers = length(fibers);

final_nodes = [(1:length(nodes))' nodes];
final_fibers = [(1:length(fibers))' fibers];

%% Crop domain into a Rectangle:
N_nodes_initial = length(final_nodes);

[final_nodes,final_fibers] = Semi_Cyl_Cut_V2(final_nodes,final_fibers,LX,LY);
% Now final_nodes and final_fibers have NaN rows where fibers and nodes
% have been cropped. Update the arrays accordingly:

[deleted_nodes,~] = find(isnan(final_nodes(:,2)));
[deleted_fibers,~] = find(isnan(final_fibers(:,2)));

kept_nodes = (1:length(final_nodes))';
kept_nodes(deleted_nodes)=[];
kept_fibers = (1:length(final_fibers))';
kept_fibers(deleted_fibers)=[];

new_node_ids = NaN(length(final_nodes),1);
new_fiber_ids = NaN(length(final_fibers),1);

id=0;
for j=kept_nodes'
    id=id+1;
    new_node_ids(j)=id;
    
end
id=0;
for j=kept_fibers'
    id=id+1;
    new_fiber_ids(j)=id;
end

%% Crop the nan rows from final nodes and fibers
final_nodes(deleted_nodes,:)=[];
N_final_nodes = length(final_nodes);
final_nodes = [(1:N_final_nodes)' final_nodes(:,2:4)];

final_fibers(deleted_fibers,:)=[];
N_final_fibers = length(final_fibers);

for j=1:N_final_fibers
    n1 = final_fibers(j,2);
    n2 = final_fibers(j,3);
    nn1 = new_node_ids(n1);
    nn2 = new_node_ids(n2);
    final_fibers(j,:) = [j nn1 nn2];
end
    


%% Create node sets for BC's
% pin all of the bottom nodes:
id = find(final_nodes(:,4)<(-(height/2)+.5));   % all nodes within .5 fiber lengths of bottom surface
boundary_nodes_2 = final_nodes(id,1);

% Apply load to opposite face nodes
id = find(final_nodes(:,4)>((height/2)-.5));
cload_nodes_2 = final_nodes(id,1);
cload_nodes_2(end)=[];

path=[cd,'\',Folder];
Orig_Dir=cd(path);
file_loc=[path,'\',cropped_network];
save(cropped_network);
cd (Orig_Dir)