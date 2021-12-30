% sanity check on nodes and fibers of cropped network:
% written by Mainak
% it removes unconnected suspending fibers

function [fibers_rev, nodes_rev] = check_raw_cropped(fibers, nodes)

nodes_int = [(1:size(nodes,1))' nodes(:,1) nodes(:,2) nodes(:,3)] ;
fibers_int = fibers ;
fibcomb = [fibers(:,1) ; fibers(:,2)] ;

count = 0 ;
for i = 1 : size(nodes_int,1)
   p = nnz(fibcomb(:,1)==nodes_int(i,1)) ;
   if p == 1 
       nodes_int(i,1:4) = [NaN NaN NaN NaN] ;
       for j = 1:size(fibers,1)
           if fibers(j,1) == i || fibers(j,2) == i 
               count = count + 1 ;
              fibers_int(j,1:2) = [NaN NaN] ;
           end
       end
   end
end

disp(['no of unconnected fibers removed = ', num2str(count)])

nodes_int2 = nodes_int ;

nodes_int2(any(isnan(nodes_int2), 2), :) = [];

renum_mat = [nodes_int2(:,1) (1:size(nodes_int2,1))'] ;

fibers_int(any(isnan(fibers_int), 2), :) = [] ;

for kk = 1:size(fibers_int,1)
    for mn = 1 : size(renum_mat,1) 
    if renum_mat(mn,1) == fibers_int(kk,1)
       fibers_int(kk,1) = renum_mat(mn,2) ;
    end 
    if renum_mat(mn,1) == fibers_int(kk,2)
       fibers_int(kk,2) = renum_mat(mn,2) ;
    end     
    end
end
    
nodes_rev = [nodes_int2(:,2) nodes_int2(:,3) nodes_int2(:,4)] ;

fibers_rev = fibers_int ;



       
