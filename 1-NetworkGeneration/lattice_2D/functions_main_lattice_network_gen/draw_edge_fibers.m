% this function can help to keep the probability of existance of fibers in the top and bottom edges of Ronceray-type
% lattice-network 
% use this function before node deletion
% written by Mainak

function [el_set_final] = draw_edge_fibers(nodes_set_final, el_set_final, edge_thkness)

side_length = max(nodes_set_final(:,3)) - min(nodes_set_final(:,3)) ;

d = edge_thkness * side_length ;   % considering some percent of the side-edge length as edge at top and bottom

j1 = find(nodes_set_final(:,3) > (side_length - d)) ;
j2 = find(nodes_set_final(:,3) < d) ;

nodes_edge1 = zeros(size(j1,1),3) ;
for i = 1 : size(j1,1)
 nodes_edge1(i,:) = nodes_set_final(j1(i),:) ;
end

nodes_edge2 = zeros(size(j2,1),3) ;
for i = 1 : size(j2,1)
 nodes_edge2(i,:) = nodes_set_final(j2(i),:) ;
end

% make the edge nodes (to apply BCs) interconnected through Delaunay
% triangulation
P1 = nodes_edge1(:,[2,3]) ;
DT1 = delaunayTriangulation(P1) ;
el_edge1 = edges(DT1) ;
el_edge1 = unique(sort(el_edge1,2),'rows') ;   % eliminate repeated fibers
el_edge1 = [el_edge1(:,1) el_edge1(:,2)] ;
for k = 1 : size(el_edge1,1)
mn = el_edge1(k,1) ;
op = el_edge1(k,2) ;
el_edge1(k,1) = nodes_edge1(mn,1) ;
el_edge1(k,2) = nodes_edge1(op,1) ;
end

P2 = nodes_edge2(:,[2,3]) ;
DT2 = delaunayTriangulation(P2) ;
el_edge2 = edges(DT2) ;
el_edge2 = unique(sort(el_edge2,2),'rows') ;   % eliminate repeated fibers
el_edge2 = [el_edge2(:,1) el_edge2(:,2)] ;
for k = 1 : size(el_edge2,1)
mn = el_edge2(k,1) ;
op = el_edge2(k,2) ;
el_edge2(k,1) = nodes_edge2(mn,1) ;
el_edge2(k,2) = nodes_edge2(op,1) ;
end

el_set  = [ el_edge1 ; el_set_final(:,2:3) ;el_edge2 ] ;
el_set = unique(sort(el_set,2),'rows') ;   % eliminate repeated fibers


index = size(el_set,1) ;

el_set_final = [(1:index)' el_set] ;



