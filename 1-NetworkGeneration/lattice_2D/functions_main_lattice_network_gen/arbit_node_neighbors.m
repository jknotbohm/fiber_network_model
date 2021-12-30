% for an arbitrary node, with node number, N, which is strictly NOT at any of the edges
% written by Mainak Sarkar

function [Neighbor_nodes] = arbit_node_neighbors(N, N1)
m = N1(N,4) ; % row number of the selected node N
n = N1(N,5) ; % column number of the selected node N
% following loop gives first neighbor node
if rem(m,2) ~= 0 
for uc = 1:size(N1, 1)
    if N1(uc,4) == (m-1) && N1(uc,5) == (n-1)
    Neighbor_nodes(1,:) =  N1(uc,:) ;
    end
end
% other neighbor nodes
for uc = 1:size(N1, 1)
    if N1(uc,4) == (m-1) && N1(uc,5) == (n)
    Neighbor_nodes(2,:) =  N1(uc,:) ;
    end
end
for uc = 1:size(N1, 1)
    if N1(uc,4) == (m) && N1(uc,5) == (n+1)
    Neighbor_nodes(3,:) =  N1(uc,:) ;
    end
end
for uc = 1:size(N1, 1)
    if N1(uc,4) == (m) && N1(uc,5) == (n-1)
    Neighbor_nodes(4,:) =  N1(uc,:) ;
    end
end
for uc = 1:size(N1, 1)
    if N1(uc,4) == (m+1) && N1(uc,5) == (n-1)
    Neighbor_nodes(5,:) =  N1(uc,:) ;
    end
end
for uc = 1:size(N1, 1)
    if N1(uc,4) == (m+1) && N1(uc,5) == (n)
    Neighbor_nodes(6,:) =  N1(uc,:) ;
    end
end
end

if rem(m,2) == 0 
for uc = 1:size(N1, 1)
    if N1(uc,4) == (m-1) && N1(uc,5) == (n)
    Neighbor_nodes(1,:) =  N1(uc,:) ;
    end
end
% other neighbor nodes
for uc = 1:size(N1, 1)
    if N1(uc,4) == (m-1) && N1(uc,5) == (n+1)
    Neighbor_nodes(2,:) =  N1(uc,:) ;
    end
end
for uc = 1:size(N1, 1)
    if N1(uc,4) == (m) && N1(uc,5) == (n+1)
    Neighbor_nodes(3,:) =  N1(uc,:) ;
    end
end
for uc = 1:size(N1, 1)
    if N1(uc,4) == (m) && N1(uc,5) == (n-1)
    Neighbor_nodes(4,:) =  N1(uc,:) ;
    end
end
for uc = 1:size(N1, 1)
    if N1(uc,4) == (m+1) && N1(uc,5) == (n)
    Neighbor_nodes(5,:) =  N1(uc,:) ;
    end
end
for uc = 1:size(N1, 1)
    if N1(uc,4) == (m+1) && N1(uc,5) == (n+1)
    Neighbor_nodes(6,:) =  N1(uc,:) ;
    end
end
end

%{
% For checking nodal coordinates
% figure
plot(N1(N,2), N1(N,3), 'X')
hold on
plot(Neighbor_nodes(:,2),Neighbor_nodes(:,3),'o')
axis equal
%}