% fiber consideration on particular node/ No of fibers linked to the
% selected arbit node, N
% written by Mainak Sarkar


function [el_set] = nodal_disorder(p, Neighbor_nodes, N1, N)      % p is the probability that the neighboring nodes are connected to the concerned node
% p = 1 implies all the six neighboring nodes are connected to the
% concerned node. p = 0.5, if 3 nodes are connected, and so on.
% p can take discrete values 3/6, 4/6, 5/6, 1.

if p == 1
    el_count = p * size(Neighbor_nodes,1) ;
el_set = [N1(N,1), Neighbor_nodes(1,1) ; 
    N1(N,1), Neighbor_nodes(2,1) ;
    N1(N,1), Neighbor_nodes(3,1) ;
    N1(N,1), Neighbor_nodes(4,1) ;
    N1(N,1), Neighbor_nodes(5,1) ;
    N1(N,1), Neighbor_nodes(6,1)] ;
end

%{
if p == 5/6 
    Neighbor_nodes_sorted = Neighbor_nodes([2,5,6],:) ;
k = randperm(size(Neighbor_nodes_sorted,1));
Neighbor_nodes_sorted = Neighbor_nodes_sorted(k(1:2),:);
el_set = [N1(N,1), Neighbor_nodes(1,1) ; 
    N1(N,1), Neighbor_nodes(3,1) ;
    N1(N,1), Neighbor_nodes(4,1) ;
    N1(N,1), Neighbor_nodes_sorted(1,1) ;
    N1(N,1), Neighbor_nodes_sorted(2,1)  ] ; 
end


if p == 4/6 
        Neighbor_nodes_sorted = Neighbor_nodes([2,5,6],:) ;
k = randperm(size(Neighbor_nodes_sorted,1));
Neighbor_nodes_sorted = Neighbor_nodes_sorted(k(1:1),:);
el_set = [N1(N,1), Neighbor_nodes(1,1) ; 
    N1(N,1), Neighbor_nodes(3,1) ;
    N1(N,1), Neighbor_nodes(4,1) ;
    N1(N,1), Neighbor_nodes_sorted(1,1) ] ;
end


if p == 3/6 
el_set = [N1(N,1), Neighbor_nodes(1,1) ; 
    N1(N,1), Neighbor_nodes(3,1) ;
    N1(N,1), Neighbor_nodes(4,1)] ;
end

%}

%{
if p == 2/6 
k = randperm(size(Neighbor_nodes,1));
Neighbor_nodes = Neighbor_nodes(k(1:2),:);
el_set = [N1(N,1), Neighbor_nodes(1,1) ; 
    N1(N,1), Neighbor_nodes(2,1)] ;
end


if p == 1/6 
k = randperm(size(Neighbor_nodes,1));
Neighbor_nodes = Neighbor_nodes(k(1:1),:);
el_set = [N1(N,1), Neighbor_nodes(1,1)] ; 
end

if p == 0 
el_set = NaN ; 
end
%}

