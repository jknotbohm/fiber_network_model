% identify and remove danging fibers, update the fiber matrix and node
% matrix
% written by Mainak Sarkar

% Goal:
% to remove the dangling fiber connected to the network only at one point

function [el_set_final, nodes_set_final] = dangling_fiber_remover(el_set_final, nodes_set_final)
nodes_to_delete = deal(NaN) ;
%
A = [el_set_final(:,2) ; el_set_final(:,3)] ;
B = nodes_set_final(:,1) ;

[GC, GR] = groupcounts(A) ;
%{
for i=1:numel(B)
    C(i) = numel(find(ismember(A, i)==1));
end
%}
D = [GR GC] ;   % D is of the format [<node number> <number of repetition>]


cn = 0 ;
for j = 1 : size(D,1) 
if D(j,2) == 1        % identify nodes that are connected to one element 
   cn = cn + 1 ;
   nodes_to_delete(cn) = D(j,1) ; % note D(j,1), nodes (identified by indices) to delete 
end 
end


TF = isnan(nodes_to_delete) ;

if TF ~= 1
%{
for k = 1 : size(nodes_to_delete,2)
nodes_set_final(nodes_to_delete(k),:) = zeros(1,3) ; 
end

nodes_set_final = nodes_set_final(any(nodes_set_final,2),:) ;
%}
    
% check all the fibers and determine the danglinng fibers, delete them
for k = 1 : size(nodes_to_delete,2)
for ws = 1:size(el_set_final, 1)
    if el_set_final(ws,2) == nodes_to_delete(k) || el_set_final(ws,3) == nodes_to_delete(k)
    el_set_final(ws,:) = zeros(1,3) ;
    end
end
end

el_set_final = el_set_final(any(el_set_final,2),:) ;

else
   el_set_final = el_set_final ;
   nodes_set_final = nodes_set_final ;
end



    



