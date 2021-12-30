% this function file corrects the numbering of the node and fiber matrix
% before feeding them into the inp generators
% written by Mainak Sarkar

function [el_set_final, nodes_set_final] = renumbering_el_nodes(el_set_final, nodes_set_final)

m = size(nodes_set_final, 1) ;

node_indx = [(1:m)' nodes_set_final(:,1)] ;

n = size(el_set_final,1) ;

for i = 1:n
    for j = 1:size(node_indx,1)
    if el_set_final(i,2) == node_indx(j,2)
        el_set_final(i,2) = node_indx(j,1) ;
    end
    end
end

for i = 1:n
    for j = 1:size(node_indx,1)
    if el_set_final(i,3) == node_indx(j,2)
        el_set_final(i,3) = node_indx(j,1) ;
    end
    end
end


nodes_set_final = [(1:m)' nodes_set_final(:,2) nodes_set_final(:,3)] ;


