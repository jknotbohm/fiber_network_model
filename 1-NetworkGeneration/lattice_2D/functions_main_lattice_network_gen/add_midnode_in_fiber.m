% this function adds a midnode to each element and revise indices accordingly
% written by Mainak Sarkar


function [el_set_final, nodes_set_final] = add_midnode_in_fiber(el_set_final, nodes_set_final)
ns = size(nodes_set_final,1) ;
ms = size(el_set_final,1) ;

for i = 1 : size(el_set_final,1)
    A_ind = el_set_final(i,2) ;
    B_ind = el_set_final(i,3) ;
    % first end point
    x1 = nodes_set_final(A_ind,2) ; 
    y1 = nodes_set_final(A_ind,3) ;
    % second end point
    x2 = nodes_set_final(B_ind,2) ; 
    y2 = nodes_set_final(B_ind,3) ;
    % midpoint of fiber
    xm = 0.5 * (x1 + x2) ;
    ym = 0.5 * (y1 + y2) ;
    % create new node and add it to the nodes_set_final matrix
    nodes_set_final((ns + i),:) = [(ns+i), xm, ym] ;  
    new_fibers_created{i} = [A_ind, (ns+i); (ns+i), B_ind] ; 
end

    
addendum_elements = cell2mat(new_fibers_created') ;
ts = size(addendum_elements,1) ;
el_set_final = [(1 : ts)' addendum_elements] ;





