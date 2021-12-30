% this function eliminates ghost nodes which are not connected to any of
% the fibers and not required to be showing up in the final analysis.
% written by Mainak Sarkar

function [final_nodes, node_index_refresher] = ghost_node_remover(final_fibers, final_nodes)
remove_node_index = NaN ;
mn = size(final_nodes,1) ;
B = final_fibers(:,2:end) ;
i = 0 ; 
for m = 1 : mn
    A = final_nodes(m,1) ;
    Lia = ismember(A,B) ;
    if Lia == 0
        i = i + 1 ;
        remove_node_index(i) = m ; 
    end
end

check = isnan(remove_node_index) ;
if check == 0
for pq = 1 : size(remove_node_index, 2)
t = remove_node_index(1, pq) ;
final_nodes(t,:) = NaN(size(final_nodes(t,:))) ;
end

final_nodes(~any(~isnan(final_nodes), 2),:)=[];

end

node_index_refresher = [(1:size(final_nodes,1))' final_nodes(:,1)] ;