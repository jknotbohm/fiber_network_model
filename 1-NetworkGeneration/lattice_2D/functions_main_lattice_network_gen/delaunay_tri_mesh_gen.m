% forms delaunay triangular lattice-based network, also Voronoi polygonal latice-based
% network
% written by Mainak Sarkar


function [nodes_set_final, el_set_final] = delaunay_tri_mesh_gen(nodes_set_final, dtvo)

if dtvo == 0 || dtvo == 1
P = nodes_set_final(:,[2,3]) ;
DT = delaunayTriangulation(P) ;
el_set_final = edges(DT) ;
el_set_final = unique(sort(el_set_final,2),'rows') ;   % eliminate repeated fibers
ind = size(el_set_final,1) ;
el_set_final = [(1:ind)' el_set_final(:,1) el_set_final(:,2)] ;

end

if dtvo == 1 
%% getting the Voronoi diagram from delaunay triangulation
[V,r] = voronoiDiagram(DT); % Compute the Voronoi vertices and regions

V = [(1:size(V,1))' V(:,1) V(:,2)] ; 
nodes_set_final = V(2:end,:) ; 


count = 0 ; 
for m = 1 : length(r)
% Obtain vertices enclosing region m
coord = V(r{m},:);
s = r{m} ;
B = any(s(:) == 1) ;
if B ~= 1
count = count + 1 ;
reg{count} = [coord(:,1), coord(:,2), coord(:,3)] ;
end
end


% getting the fiber matrix
count2 = 0 ; 
for mn = 1 : count
    lattice_el = reg{mn} ;
    no_of_sides = size(lattice_el,1) ;
    for uz = 1:no_of_sides
        if uz ~= no_of_sides
        count2 = count2 + 1 ;
       fibers{count2} = [lattice_el(uz,1), lattice_el(uz+1,1) ] ;
        else
            count2 = count2 + 1 ;
            fibers{count2} = [lattice_el(uz,1), lattice_el(1,1) ] ;
        end
    end
end

fibers = cell2mat(fibers') ;

%{
% eliminate repeated fibers
for kk = 1:size(fibers,1) 
    for ll = 1:size(fibers,1)
    if fibers(kk,1) == fibers(ll,2)
        if fibers(ll,1) == fibers(kk,2)
           fibers(kk,:) = zeros(1, size(fibers,2)) ;
        end
    end
    end
end
%}  
      
% fibers = fibers(any(fibers,2),:) ;

% eliminate repeated fibers (more efficient)
fibers = unique(sort(fibers,2),'rows');

el_set_final = [(1:size(fibers,1))' fibers(:,1) fibers(:,2)] ; 


end


