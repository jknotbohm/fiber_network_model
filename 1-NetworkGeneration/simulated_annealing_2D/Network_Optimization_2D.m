
function [ nodes, fibers, fiberenergy, fiberlengths, total_energy, N1,N2,N_interior ] = Network_Optimization_2D( fraction_to_try_swap,N, nodes, fibers, N_anneal, lx,ly,l_fiber,fiberlengths,fiberenergy,N1,N2,N_boundary_nodes, stepsize,swap_skip_energy )
% Network_Optimization: use simulated annealing to iterate initial network

% Parameters
N_interior = [N_boundary_nodes+1:N]';   % The boundary nodes are at the top of nodes, we don't want to move those
N_int = length(N_interior);
%% ---- Network Optimization: -----

% Initialize counters:
anneal_id = 0;
swap_id=0;
n_swaps_accepted = 0;
n_swaps_rejected = 0;
n_accepted = 0;
n_prob=0;

%Initialize node_fiber_matrix
node_fiber_mat = NaN(length(nodes),10);

for i=1:length(nodes)
    [nlf,~]=find(fibers==i);
    node_fiber_mat(i,1:length(nlf))=nlf';
end

for m=1:N_anneal
    fork = rand;
    if fork > fraction_to_try_swap
        % Do node displacement, just on interior nodes.
        for j=N_interior'
            anneal_id = anneal_id+1;
            % Try displacing the jth node by a random 3D spatial step:
            newx = nodes(j,1) + stepsize*(2*(-.5 + rand))*(((N_anneal*N_int)-n_accepted)/(N_anneal*N_int));
            newy = nodes(j,2) + (ly/lx)*stepsize*(2*(-.5 + rand))*(((N_anneal*N_int)-n_accepted)/(N_anneal*N_int));
            
            % Find the fibers that this would affect:
            n1f = node_fiber_mat(j,~isnan(node_fiber_mat(j,:)));
            n1f=n1f';
            othernodes = zeros(length(n1f),1);
            for k=1:length(n1f)
                node1 = fibers(n1f(k),1);
                node2 = fibers(n1f(k),2);
                if node1 == j
                    othernodes(k) = node2;
                else
                    othernodes(k) = node1;
                end
            end
            
            nfl = zeros(length(n1f),1); % new fiber lengths
            for k=1:length(nfl)
                nfl(k) =  ((nodes(othernodes(k),1)-newx)^2 + (nodes(othernodes(k),2)-newy)^2)^.5;
            end
            
            nfe = (nfl-l_fiber).^2;
            nfe_sum = sum(nfe);
            % compare this to the sum of the corresponding fiberenergy
            % values.
            
            current_fes = fiberenergy(n1f);
            current_fes_sum = sum(current_fes);
            
            if nfe_sum < current_fes_sum
                nodes(j,:) = [newx,newy];      % update relevant arrays
                
                % For the indices (cni), either N1 or N2 corresponds to the new
                % node position. Now that nodes has been updated we can just
                % update both N1 and N2 of cni even though only one of each
                % actually has changed.
                
                N1(n1f,:) = nodes(fibers(n1f,1),:);
                N2(n1f,:) = nodes(fibers(n1f,2),:);
                
                fiberlengths(n1f) = ((N1(n1f,1)-N2(n1f,1)).^2 + (N1(n1f,2)-N2(n1f,2)).^2 ).^.5;
                fiberenergy(n1f) = (fiberlengths(n1f) - l_fiber).^2;
                n_accepted = n_accepted + 1;
            elseif nfe_sum < 1.2*current_fes_sum && rand(1)>.95   % Let there be a small percent of slightly positive energy transitions to help get out of local minima
                nodes(j,:) = [newx,newy];      % update relevant arrays
                
                % For the indices (cni), either N1 or N2 corresponds to the new
                % node position. Now that nodes has been updated we can just
                % update both N1 and N2 of cni even though only one of each
                % actually has changed.
                
                N1(n1f,:) = nodes(fibers(n1f,1),:);
                N2(n1f,:) = nodes(fibers(n1f,2),:);
                
                fiberlengths(n1f) = ((N1(n1f,1)-N2(n1f,1)).^2 + (N1(n1f,2)-N2(n1f,2)).^2).^.5;
                fiberenergy(n1f) = (fiberlengths(n1f) - l_fiber).^2;
                n_accepted = n_accepted + 1;
            end
        end
        
    else
        % Do node swapping: this procedure finds two fibers, fiber1 and
        % fiber2, whose nodes are called [node1 node3] and [node2 node4]
        % and swaps it, so fiber1 becomes [node1 node2] and fiber2 becomes
        % [node3 node4]
        for j=1:N

            node1 = j;
            % Find node 1's fibers:
%             [n1f,~] = find(fibers==node1);  % fibers that use node1
            n1f=node_fiber_mat(node1,~isnan(node_fiber_mat(node1,:)));
            n1f=n1f';
            
            current_nodal_energy = fiberenergy(n1f);
            cne_sum = sum(current_nodal_energy);
            
            if cne_sum < swap_skip_energy
                continue            % skip this node if its energy is already pretty low.
            end
            
            anneal_id = anneal_id+1;
            swap_id=swap_id+1;
            
            % Find the nodes of the fibers connected to the jth node          
            node1fibers_nodes=zeros(length(n1f),1);
            for h=1:length(node1fibers_nodes)
                n1 = fibers(n1f(h),1);
                n2 = fibers(n1f(h),2);
                if n1 == node1
                    node1fibers_nodes(h) = n2;
                elseif n2 == node1
                    node1fibers_nodes(h) = n1;
                end
            end
            
            % Choose which fiber to delete:            
            [~,id]=max(current_nodal_energy);
            
            nftd_1 = node1fibers_nodes(id);  % pick the worst fiber of node j
            ftd_1 = n1f(id);

            distances_1 = ((nodes(node1,1) - nodes(:,1)).^2 + (nodes(node1,2) - nodes(:,2)).^2 ).^.5;
            neighbors_1 = find(distances_1<2*l_fiber & distances_1>.5*l_fiber);   % node indices
            
            % Throw out the nodes already connected to the jth node
            for k=1:length(node1fibers_nodes)
                cfsin = find(neighbors_1==node1fibers_nodes(k));
                neighbors_1(cfsin)=[];
            end
            
            % Now we also need to throw out any neighbors that are
            % connected to nftd_1:

            nftd_1_fibers = node_fiber_mat(nftd_1,~isnan(node_fiber_mat(nftd_1,:)));
            nftd_1_fibers=nftd_1_fibers';
            nftd_1_nodes=zeros(length(nftd_1_fibers),1);
            for h=1:length(nftd_1_fibers)
                n1 = fibers(nftd_1_fibers(h),1);
                n2 = fibers(nftd_1_fibers(h),2);
                if n1 == nftd_1
                    nftd_1_nodes(h) = n2;
                elseif n2 == nftd_1
                    nftd_1_nodes(h) = n1;
                else
                    disp('nftd_1 node search is busted');
                end
            end
            
            for k=1:length(nftd_1_nodes)
                cfsin = find(neighbors_1==nftd_1_nodes(k));
                neighbors_1(cfsin)=[];
            end
            
            % Now any remaining nodes in neighbors_1 should be viable swap
            % candidates
            
            % If there's none left continue on to next iteration:
            
            if length(neighbors_1)==0
                continue
            end
            
            node3=nftd_1;
            
            % Choose a node from neighbors_1, find its worst fiber:
            wc=0;
            while wc < length(neighbors_1)
                wc=wc+1;
                n2try = neighbors_1(wc);
                node2=n2try;
                % NEW CODE 
                n2try_fibers = node_fiber_mat(n2try,~isnan(node_fiber_mat(n2try,:)));
                n2try_fibers = n2try_fibers';
                if isempty(n2try_fibers)
                    continue
                end

                current_nodal_energy = fiberenergy(n2try_fibers);
                [~,id]=max(current_nodal_energy);
                fiber2 = n2try_fibers(id);
                
                % Find the other node of fiber2
                n1=fibers(fiber2,1);    % dummies
                n2=fibers(fiber2,2);
                if n1==node2
                    node4=n2;
                else
                    node4=n1;
                end
                
                % check that this node isn't connected to node 1 or node 3:

                node4check = find(node1fibers_nodes==node4);
                if length(node4check)>0;
                    n_prob=n_prob+1;
                    continue
                end
                
                node4check = find(nftd_1_nodes==node4);
                if length(node4check)>0;
                    n_prob=n_prob+1;
                    continue
                end                
                
                
                % old code that went here pasted on bottom of script
                
                fiber1 = ftd_1;

                % Find the two new nodes associated with this:
                
                % currently fiber1 = [node1, node3] and fiber2 = [node2, node4]
                
                node1 = fibers(fiber1,1);
                node2 = fibers(fiber1,2);
                node3 = fibers(fiber2,1);
                node4 = fibers(fiber2,2);
                
                if node1==node4;
                    continue
                elseif node3==node2;
                    continue
                end
                
%                 if n1 == node1
%                     node3 = n2;
%                 else
%                     node3 = n1;
%                 end
%                 
%                 if n3 == node2
%                     node4 = n4;
%                 else
%                     node4 = n3;
%                 end
                
                % Now swap it so that fiber1 = [node1,node4] and fiber2 =
                % [node2,node3]
                
                new_fiber_1 = [node1,node4];
                new_fiber_2 = [node2,node3];
                
                
                % compare this to the actual fiber1 and fiber2:
                
                newfiber1length = ((nodes(node1,1)-nodes(node4,1))^2 + (nodes(node1,2)-nodes(node4,2))^2 )^.5;
                newfiber2length = ((nodes(node2,1)-nodes(node3,1))^2 + (nodes(node2,2)-nodes(node3,2))^2 )^.5;
                
                newfibersenergy = (newfiber1length - l_fiber)^2 + (newfiber2length - l_fiber)^2;
                oldtwofibers = [fiber1;fiber2];
                oldfibersenergy = sum(fiberenergy(oldtwofibers));
                
                % Now check that we didn't create a duplicate fiber in the node
                % swapping process:
                %             newfibers = fibers;
                %             newfibers(fiber1,:) = new_fiber_1;
                %             newfibers(fiber2,:) = new_fiber_2;
                %
                %             N_new_fibers=length(newfibers);
                %             setOrder='stable';
                %             fdc = [newfibers(:,2) newfibers(:,1)];
                %             fiber_duplicates = [newfibers; fdc];
                %             [~,fiber_duplicate_check,~] = unique(fiber_duplicates,'rows',setOrder);
                %             N_fiber_duplicates = 2*N_new_fibers - length(fiber_duplicate_check);
                
                
                if newfibersenergy < oldfibersenergy %&& N_fiber_duplicates == 0
                    
                    %remove fib in node_fiber_mat
                    for n=1:2
                        node_rem = fibers(fiber1,n);
                        fib_rem_ind=find( node_fiber_mat(node_rem,:)== fiber1 );
                        node_fiber_mat(node_rem,fib_rem_ind(1))=NaN;
                    end
                    
                    for n=1:2
                        node_rem = fibers(fiber2,n);
                        fib_rem_ind=find( node_fiber_mat(node_rem,:)== fiber2 );
                        node_fiber_mat(node_rem,fib_rem_ind(1))=NaN;
                    end
                    
                    %Change actual Fiber location
                    fibers(fiber1,:) = new_fiber_1;
                    fibers(fiber2,:) = new_fiber_2;
                    
                    %add fib to node_fiber_mat                    
                    for n=1:2
                        node_add = new_fiber_1(n);
                        fib_add_ind=find( isnan(node_fiber_mat(node_add,:)));
                        node_fiber_mat(node_add,fib_add_ind(1))= fiber1;
                    end
                    
                    for n=1:2
                        node_add = new_fiber_2(n);
                        fib_add_ind=find( isnan(node_fiber_mat(node_add,:)));
                        node_fiber_mat(node_add,fib_add_ind(1))= fiber2;
                    end
                    
                    N1(fiber1,:) = nodes(node1,:);
                    N2(fiber1,:) = nodes(node4,:);
                    N1(fiber2,:) = nodes(node2,:);
                    N2(fiber2,:) = nodes(node3,:);
                    
                    fiberlengths(n1f) = ((N1(n1f,1)-N2(n1f,1)).^2 + (N1(n1f,2)-N2(n1f,2)).^2 ).^.5;
                    fiberenergy(n1f) = (fiberlengths(n1f) - l_fiber).^2;
                    
                    n_accepted = n_accepted+1;
                    n_swaps_accepted = n_swaps_accepted+1;
                    break;  % exit the while loop
                else
                    n_swaps_rejected = n_swaps_rejected+1;
                    %wc=wc+1;
                end
            end
        end
    end
    
end

percent_accepted_iterations = 100*n_accepted/anneal_id;
str=['Percent accepted iterations for that optimization run = ' num2str(percent_accepted_iterations)];
disp(str);
%percent_accepted_swaps = 100*n_swaps_accepted/swap_id;

total_energy = sum(fiberenergy);
str=['Total Network Length Energy = ' num2str(total_energy)];
disp(str);
avg_fiber_length = mean(fiberlengths);
median_fiber_length = median(fiberlengths);
str = ['mean fiber length = ' num2str(avg_fiber_length) ' median length = ' num2str(median_fiber_length)];
disp(str);
N_fibers=length(fibers);
setOrder='stable';
fdc = [fibers(:,2) fibers(:,1)];
fiber_duplicates_final = [fibers; fdc];
[~,fiber_duplicate_check_final,~] = unique(fiber_duplicates_final,'rows',setOrder);
final_fiber_check = 2*N_fibers - length(fiber_duplicate_check_final);

str = ['The number of duplicate fibers is ' num2str(final_fiber_check)];
disp(str);

end
