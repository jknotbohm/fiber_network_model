function [ nodes, fibers, fiberenergy, total_energy, N1, N2, N_boundary_nodes,fiberlengths,valencycheck ] = Initial_Network_Generation_2D(rho,lx,ly,l_fiber,n_neighbors_to_try)
% Initial_Network_Generation: Create an initial fiber network configuration
% that satisfies the desired valency distribution, but not necessarily
% other network properties.

% Note: n_neighbors_to_try is the number of the closest nodes to each node
% to attempt to connect to in trying to achieve the desired valency for
% each node. If n_neighbors_to_try is too low, then there will be a large
% amount of nodes with the wrong valency since they won't be able to find
% connecting nodes that aren't already saturated. If n_neighbors_to_try is
% too large, then the initial network will have fibers that are excessively
% long, and it will take longer for optimization to converge.

N = round(rho*lx*ly);    % Number of nodes

NBx = round(.5*(N/lx));
NBy = round(.5*(N/ly));  % These are the number of nodes that will be placed on the x, y, and z face
% NBz = round(.5*(N/lz));  % boundaries. These are fixed while the interior node positions are optimized

%% ---- INITIAL NETWORK GENERATION-----

% Boundary nodes:

N_boundary_nodes = 2*(NBx+NBy);

nodes_boundary = zeros(N_boundary_nodes,2);

% X boundary:

xbnd1_x = (-lx/2)*ones(NBx,1);
xbnd1_y = (-ly/2) + ly*rand(NBx,1); 

xbnd2_x = (lx/2)*ones(NBx,1);
xbnd2_y = (-ly/2) + ly*rand(NBx,1); 

nodes_boundary(1:2*NBx,:) = [xbnd1_x xbnd1_y; xbnd2_x xbnd2_y];

ybnd1_x = (-lx/2) + lx*rand(NBy,1);
ybnd1_y = (-ly/2)*ones(NBy,1);

ybnd2_x = (-lx/2) + lx*rand(NBy,1);
ybnd2_y = (ly/2)*ones(NBy,1);

nodes_boundary(2*NBx+1:2*NBx+2*NBy,:) = [ybnd1_x ybnd1_y ; ybnd2_x ybnd2_y ];


%% Preallocate node array:
nodes = zeros(N,2); %(nodeindex, X, Y, Z)

nodes(1:N_boundary_nodes,:) = nodes_boundary;

% The rest of the nodes should be uniformly randomly distributed inside the
% domain.
nodes(N_boundary_nodes+1:N,:) = ...
    rand(length(nodes)-N_boundary_nodes,2).*[lx ly] - 0.5*[lx ly];

% Each node needs to be assigned a valency of either 3 or 4. Build a valency array:

% In the future there should be some sort of valency probability
% distribution and we can select each nodal valency by using rand and
% rounding to the nearest valency value.

valency = rand(N,1);
valency(valency < .8)=3;
valency(valency < .9)=4;
valency(valency < .97)=5;
valency(valency < 1)=6;


% Now, connect nodes with fibers such that the assigned valency
% for each node is satisfied. The initial network needs to be somewhat
% realistic so only connect nodes to other close neighbors, don't just
% randomly connect nodes together throughout the domain.

% Preallocate the fibers array (don't know how big it will be so make an
% oversized zeros array and we'll crop it later)

fibers = zeros(4*N,2);      % each row of fibers defines a fiber with two node numbers
current_valency_matrix=zeros(N,2); current_valency_matrix(:,1)=1:N;
fiberid = 1;
n_wll=0;
for k=1:N
    %Check current valency of the k'th node by seeing how many times its node
    %number appears in a row of the fibers array:
    currentvalency=current_valency_matrix(k,2);
    
    node_loc = nodes(k,:);
    neighbors=0; node_dist=5;
    while length(neighbors) < n_neighbors_to_try+1
        node_dist=node_dist+5;
        x_lim = find((nodes(:,1)>node_loc(1)-node_dist) & (nodes(:,1)<node_loc(1)+node_dist));
        y_lim = find((nodes(:,1)>node_loc(1)-node_dist) & (nodes(:,1)<node_loc(1)+node_dist));

        neighbors=intersect(x_lim,y_lim);
        
    end
    
    distances = ((nodes(k,1) - nodes(neighbors,1)).^2 + (nodes(k,2) - nodes(neighbors,2)).^2 ).^.5;
    [~, dis] = sort(distances);
    
    neighbors = neighbors(dis(2:n_neighbors_to_try+1));      % 1st dis is the current node
    %while currentvalency < valency(k)
    % Throw out any neighbors that are already connected to the k'th node
    
    for m=1:(valency(k)-currentvalency)     % add however many fibers to this node as are still needed.
        fibers(fiberid,1) = k;
        current_valency_matrix(k,2) = current_valency_matrix(k,2)+1;
        eid = randi(n_neighbors_to_try);
        endpoint = neighbors(eid);          % pick a random node from neighbors to use as endpoint
        
        %Now check the valency of this node as well:
        currentneighborvalency=current_valency_matrix(endpoint,2);
        whilelooplimit = 10*n_neighbors_to_try;   % avoid infinite loop
        wll=0;
        while currentneighborvalency >= valency(endpoint)
            if wll > whilelooplimit
                n_wll = n_wll+1;
                break               % can't always satisfy valency.
            end
            %eid=eid+1;
            eid = randi(n_neighbors_to_try);
            endpoint = neighbors(eid);                       % Can't connect to a node already at its valency
            currentneighborvalency=current_valency_matrix(endpoint,2);
            wll=wll+1;
        end
        fibers(fiberid,2) = endpoint;
        current_valency_matrix(endpoint,2) = current_valency_matrix(endpoint,2)+1;
        fiberid = fiberid+1;
        
    end
    
    
end

%crop fibers:

fibers(~any(fibers,2),:) = [];

% Remove any fibers that have the same node as both endpoints (zero
% length). Cannot prescribe that zero length fibers can't be created as
% they need to be able to be created for the initial network generation
% scheme to work. They are undesirable but as long as nzerolength/N is low
% it's not a big deal overall.

zlc = fibers(:,1)-fibers(:,2);
zlcheck = find(zlc==0);
fibers(zlcheck,:)=[];

% Check for literal duplicates in fibers:

setOrder = 'stable';
[~,fiber_duplicate_check,~] = unique(fibers,'rows',setOrder);

% Only keep the unique rows of fibers:

fibers = fibers(fiber_duplicate_check,:);
N_fibers = length(fibers);

% Now check for "mirrored" duplicates, e.g. two rows of fibers: [1 534; 534 1]

fmc = [fibers(:,2) fibers(:,1)];
fibers_plus_mirrored = [fibers; fmc];
setOrder='stable';
[~,fiber_mirror_check,~] = unique(fibers_plus_mirrored,'rows',setOrder);

% fibers_plus_mirrored(1:N_fibers) is the same as fibers, b/c there are no
% duplicates in that section of the array.

mirrors_deleted = fibers_plus_mirrored(fiber_mirror_check,:);
crop_off_fibers = mirrors_deleted(N_fibers+1:end,:);
fibers_with_mirrors_deleted = [crop_off_fibers(:,2) crop_off_fibers(:,1)];

fibers = fibers_with_mirrors_deleted;

%Check what the average valency is:

valencycheck = zeros(N,1);

for k=1:N
    [vcheck,~] = find(fibers==k);
    valencycheck(k) = length(vcheck);
end

val_0 = find(valencycheck==0);
str=['N valency=0 is ' num2str(length(val_0))];
disp(str);
val_1 = find(valencycheck==1);
str=['N valency=1 is ' num2str(length(val_1))];
disp(str);
val_2 = find(valencycheck==2);
str=['N valency=2 is ' num2str(length(val_2))];
disp(str);

add_fibers_0_1 = zeros(length(val_0),2);
add_fibers_0_2 = add_fibers_0_1;
add_fibers_0_3 = add_fibers_0_1;
add_fibers_1_1 = zeros(length(val_1),2);
add_fibers_1_2 = add_fibers_1_1;
add_fibers_3 = zeros(length(val_2),2);

%fibers = [fibers;add_fibers];    % preallocate

for k=1:length(val_0)
    % Need to add two more fibers to each of the val_1 nodes
    distances = ((nodes(val_0(k),1) - nodes(:,1)).^2 + (nodes(val_0(k),2) - nodes(:,2)).^2 + (nodes(val_0(k),3) - nodes(:,3)).^2).^.5;
    [~, dis] = sort(distances);
    neighbor = dis(2:4); % going to add one fiber between each of the two closest neighbors.
    add_fibers_0_1(k,:)=[val_0(k) neighbor(1)];
    add_fibers_0_2(k,:)=[val_0(k) neighbor(2)];
    add_fibers_0_3(k,:)=[val_0(k) neighbor(3)];
end

for k=1:length(val_1)
    % Need to add two more fibers to each of the val_1 nodes
    distances = ((nodes(val_1(k),1) - nodes(:,1)).^2 + (nodes(val_1(k),2) - nodes(:,2)).^2 ).^.5;
    [~, dis] = sort(distances);
    neighbor = dis(2:3); % going to add one fiber between each of the two closest neighbors.
    add_fibers_1_1(k,:)=[val_1(k) neighbor(1)];
    add_fibers_1_2(k,:)=[val_1(k) neighbor(2)];
end

for k=1:length(val_2)
    % Need to add two more fibers to each of the val_1 nodes
    distances = ((nodes(val_2(k),1) - nodes(:,1)).^2 + (nodes(val_2(k),2) - nodes(:,2)).^2 ).^.5;
    [~, dis] = sort(distances);
    neighbor = dis(2); % going to add one fiber between each of the two closest neighbors.
    add_fibers_3(k,:)=[val_2(k) neighbor(1)];
end

fibers = [fibers;add_fibers_0_1;add_fibers_0_2;add_fibers_0_3;add_fibers_1_1;add_fibers_1_2;add_fibers_3];

for k=1:N
    [vcheck,~] = find(fibers==k);
    valencycheck(k) = length(vcheck);
end

check = valency - valencycheck;
wrong = find(check);

n_wrong_valency = length(wrong);
avgvalency = mean(valencycheck);

str = ['The average valency is ' num2str(avgvalency)];
disp(str);

wrongvalencyfraction=n_wrong_valency/N;       % Want this to be low
str=['The percent of nodes with the wrong valency is ' num2str(100*wrongvalencyfraction)];
disp(str);

% Remove any duplicate fibers:

figure
histogram(valencycheck);
title('node valency distribution');

% Set up arrays with the fiber endpoint coordinates.
N_fibers = length(fibers);
N1 = zeros(N_fibers,2);
N2 = N1;

for m=1:N_fibers
    
    N1(m,:) = nodes(fibers(m,1),:);
    N2(m,:) = nodes(fibers(m,2),:);
    
end

fiberlengths = ((N1(:,1)-N2(:,1)).^2 + (N1(:,2)-N2(:,2)).^2 ).^.5;

% Find the total_energy of the initial network:
fiberenergy = (fiberlengths - l_fiber).^2;
total_energy = sum(fiberenergy);

end

