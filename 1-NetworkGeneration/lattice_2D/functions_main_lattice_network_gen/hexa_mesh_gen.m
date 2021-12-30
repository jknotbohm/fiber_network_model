% generating hexagonal mesh from the random triangular mesh
% written by Mainak Sarkar

function [el_set_final] = hexa_mesh_gen(N1, Nx, Ny, el_set_final, nodes_set_final)

%%
for m = 3 : 2 : (Ny - 2)
for n = 1 : 3 : (Nx - 0)
    
for t = 1 : size(N1,1)
            if N1(t,4) == m && N1(t,5) == n 
                concerned_node(1,:) = nodes_set_final(t,:) ;
            end
end
% find the number of fibers attached to that concerned node
cnt1 = 0 ;
for v = 1 : size(el_set_final,1)
    if el_set_final(v,2) == concerned_node(1,1) || el_set_final(v,3) == concerned_node(1,1)
        cnt1 = cnt1+1 ; % counts fibers connected to the concerned node
        fiber_index_current_node(cnt1) = v ;
    end
end

for vv = 1 : size(el_set_final,1)
for uu = 1:size(fiber_index_current_node,2)
    if fiber_index_current_node(uu) == el_set_final(vv,1)
        el_set_final(vv,:) = zeros(1,3) ; 
    end
end
end

end
end


%%
for m = 4 : 2 : (Ny - 2)    % Ny - 3
for n = 2 : 3 : (Nx - 0)    % Nx - 4
    
for t = 1 : size(N1,1)
            if N1(t,4) == m && N1(t,5) == n 
                concerned_node(1,:) = nodes_set_final(t,:) ;
            end
end
% find the number of fibers attached to that concerned node
cnt1 = 0 ;
for v = 1 : size(el_set_final,1)
    if el_set_final(v,2) == concerned_node(1,1) || el_set_final(v,3) == concerned_node(1,1)
        cnt1 = cnt1+1 ; % counts fibers connected to the concerned node
        fiber_index_current_node(cnt1) = v ;
    end
end

for vv = 1 : size(el_set_final,1)
for uu = 1:size(fiber_index_current_node,2)
    if fiber_index_current_node(uu) == el_set_final(vv,1)
        el_set_final(vv,:) = zeros(1,3) ; 
    end
end
end

end
end

%%
el_set_final = el_set_final(any(el_set_final,2),:) ;
