% following function will find the average connectivity of the present
% network (will leave the edge portion)
% written by Mainak

function [average_nodal_connectivity] = find_avg_node_connectivity(el_set_final, nodes_set_final, edge_thkness)

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


nodes_rejection = [nodes_edge1(:,1) ; nodes_edge2(:,1)] ;


for k = 1 : size(el_set_final,1)
    if el_set_final(k,1) ~= 0 
       u1 = ismember(el_set_final(k,2),nodes_rejection) ;
       u2 = ismember(el_set_final(k,3),nodes_rejection) ;
if u1 == 1 && u2 == 1
    el_set_final(k,:) = zeros(1,3) ;
end
    end
end

el_set_final = el_set_final(any(el_set_final,2),:) ;


B = [ el_set_final(:,2) ; el_set_final(:,3) ] ;

uv = unique(B);
n  = histcounts(B,uv);

average_nodal_connectivity = mean(n) ;






