% this will randomly delete some fibers as per probability 'prob' of
% existance of a fiber. 
% written by Mainak Sarkar

function [el_set_final] = random_fiber_eliminator_adv(pr, el_set_final, nodes_set_final, ds)

side_length = max(nodes_set_final(:,3)) - min(nodes_set_final(:,3)) ;

d = ds * side_length ;   % considering some percent of the side-edge length as edge at top and bottom

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
if u1 == 0 && u2 == 0
% check whether each fiber is existing or not
q = rand;
if q > pr
y = 0 ; % zero means fiber needs to be deleted
else
    y = 1 ; % retain the fiber 
end
       
if y == 0 
    % delete the fiber concerned
    el_set_final(k,:) = zeros(1,3) ;
end
end
    end
end


el_set_final = el_set_final(any(el_set_final,2),:) ;

    

