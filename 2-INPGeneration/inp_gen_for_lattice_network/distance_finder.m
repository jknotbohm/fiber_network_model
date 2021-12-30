% This script calculates average fiber length in the network
% written by Mainak Sarkar


function [l_fiber] = distance_finder(nodes_set_final, el_set_final)

for w = 1 : size(el_set_final,1)
    m = el_set_final(w,2) ; n = el_set_final(w,3) ; 
X = [nodes_set_final(m,2:3) ; nodes_set_final(n,2:3)];
d(w) = pdist(X,'euclidean') ;
end

l_fiber = mean(d) ;