%% Branch Optimization

function [ nodes, fibers, nodal_branching_energy, total_branching_energy_init, total_branching_energy_final ] = Branch_Optimization( N_branching_optimize, nodes, fibers, N );

% This script will align all the fibers connected at a node such that the
% two straightest fibers are straightened, but all other fibers are aligned
% towards one of the two straight fibers. This should better represent the
% physical gels' branching.

stepsize_mag = .15;

othernodes = cell(1,N);
valency = zeros(N,1);
n2s = zeros(N,1);

for j=1:N
    node1=j;
    [cf,~] = find(fibers==node1);
    n_cf = length(cf);

    valency(j) = n_cf;
    
    othernodes_jth = zeros(n_cf,1);
    
    for k=1:n_cf
        n1 = fibers(cf(k),1);
        n2 = fibers(cf(k),2);
        if n1 == node1
            othernodes_jth(k) = n2;
        else
            othernodes_jth(k) = n1;
        end
    end
    
    othernodes{1,node1} = othernodes_jth;
end

% othernodes now has an array for each node that lists the other nodes that
% connect to it via fibers.

%% Optimize the branching energy of each node

% Formulate the initial branching energies:

nodal_branching_energy = zeros(N,1);
% store the two straightest fibers at the node:
sf1 = zeros(N,1);
sf2 = sf1;

% store some things that will be getting continuously updated:

of_stored = cell(N,2);
unit_vecs = cell(N,1);
dot_prods = cell(N,1);
SFs = cell(N,1);
fiber_to_align_toward = zeros(N,1);

for j=1:N
    node1=j;
    othernodes_jth = othernodes{1,node1};
    n_cf = length(othernodes_jth);
    unit_vecs{j} = zeros(n_cf,3);
    cf_lengths = zeros(n_cf,1);
    for k=1:n_cf
        cf_lengths(k) = ((nodes(othernodes_jth(k),1)-nodes(node1,1))^2 + (nodes(othernodes_jth(k),2)-nodes(node1,2))^2 + (nodes(othernodes_jth(k),3)-nodes(node1,3))^2)^.5;
        unit_vecs{j}(k,1) = (nodes(othernodes_jth(k),1) - nodes(node1,1))/cf_lengths(k);
        unit_vecs{j}(k,2) = (nodes(othernodes_jth(k),2) - nodes(node1,2))/cf_lengths(k);
        unit_vecs{j}(k,3) = (nodes(othernodes_jth(k),3) - nodes(node1,3))/cf_lengths(k);
    end
    
    % Now find the dot products between all of the unit vectors of the
    % fibers that are coming out of this node
    
    N_combinations = nchoosek(n_cf,2);
    dot_prods{j} = zeros(N_combinations,1);
    combos = combnk(1:n_cf,2);
    
    for k=1:N_combinations
        dot_prods{j}(k) = dot(unit_vecs{j}(combos(k,1),:),unit_vecs{j}(combos(k,2),:));
    end
    
    [mindot,id]=min(dot_prods{j});
    SFs{j} = combos(id,:);
    sf1 = SFs{j}(1);   sf2 = SFs{j}(2);
    
    other_fibers = combos;
    other_fibers(id,:)=[];     % in case there are only 2 fibers at a node, maybe we don't want to put a []
                                % here because that would mess up the
                                % indexing later?
    ofs = [1:length(othernodes_jth)]';
    id2 = [sf1;sf2];
    ofs(id2)=[];
    of_stored{j} = ofs;
    %of_stored{j} = other_fibers;
    
    term_1 = 1+mindot;
    
    % sf1 and sf2 are the two straightest fibers among the ones at this
    % node. These fibers should be made as straight as possible while all
    % the other fibers need to align towards either sf1 or sf2.
    
    % Dot sf1 with all the others, and then dot sf2 with all the others,
    % and see which one has closer alignment
    
    other_unit_vecs = unit_vecs{j};
    other_unit_vecs([sf1;sf2],:)=[];
    %other_unit_vecs(sf2,:)=[];
    sf1_uv = unit_vecs{j}(sf1,:);
    sf2_uv = unit_vecs{j}(sf2,:);
    
    dot_prods_sf1 = zeros(size(other_unit_vecs,1),1);
    dot_prods_sf2=dot_prods_sf1;
    
    for k=1:size(other_unit_vecs,1);
        dot_prods_sf1(k) = dot(sf1_uv,other_unit_vecs(k,:));
        dot_prods_sf2(k) = dot(sf2_uv,other_unit_vecs(k,:));
    end
    
    align_sf1 = sum(dot_prods_sf1);
    align_sf2 = sum(dot_prods_sf2);
    
    if align_sf1 > align_sf2
        term_2 = size(other_unit_vecs,1)-align_sf1;
        fiber_to_align_toward(j,1)=1;
    else
        term_2 = size(other_unit_vecs,1)-align_sf2;
        fiber_to_align_toward(j,2)=1;
    end
    
    nodal_branching_energy(j,1)= 8*term_1 + term_2;
end
figure
histogram(nodal_branching_energy)
title('Nodal Branching Energy Before Branching Optimization')

total_branching_energy_init = sum(nodal_branching_energy);

% Now we want to try looping through all the nodes and randomly displacing
% them to minimize the nodal_branching_energy.

N_accepted_BO = 0;
N_accepted_wrt_m = zeros(N_branching_optimize,1);
total_BE_intermed=N_accepted_wrt_m;
opt_id = 0;
n_baa=0;
H=waitbar(0,'Branch optimizing...');
for m=1:N_branching_optimize
    stepsize = stepsize_mag;%*((N_branching_optimize-m+1)/N_branching_optimize);
    waitbar(m/N_branching_optimize);
    for j=1:N
        opt_id=opt_id+1;
        node1=j;
        oldx = nodes(node1,1);
        oldy = nodes(node1,2);
        oldz = nodes(node1,3);
        
        n_cf = length(othernodes_jth);        

%         if n_cf==1      % in case a valency=1 node somehow made it through
%             break
%         end
        
        endpt = (unit_vecs{j}(SFs{j}(1),:) + unit_vecs{j}(SFs{j}(2),:));
        l_endpt = sqrt((endpt(1))^2 + (endpt(2))^2 + (endpt(3))^2);
        
        dir_to_displace = endpt/l_endpt;      % now a unit vector
        
        newx = oldx + stepsize*dir_to_displace(1) + (stepsize/50)*(-.5+rand(1));
        newy = oldy + stepsize*dir_to_displace(2) + (stepsize/50)*(-.5+rand(1));
        newz = oldz + stepsize*dir_to_displace(3) + (stepsize/50)*(-.5+rand(1));
        
        % Try displacing the othernodes connected to it randomly also
        
        othernodes_jth = othernodes{node1};
        n_cf = length(othernodes_jth);        

        cnis = [node1;othernodes_jth];
        other_old_coords = nodes(othernodes_jth,:);
        other_spatial_steps = zeros(n_cf,3);
        
        for k=1:(n_cf-2)
            
            endpt = unit_vecs{j}(SFs{j}(find(fiber_to_align_toward(j,:))),:) - unit_vecs{j}(of_stored{j}(k),:);
            
            %endpt = unit_vecs{j}(SFs{j}(find(fiber_to_align_toward(j,:))),:) - unit_vecs{j}(of_stored{j}(k),:);
            %dir_to_displace = endpt;      % let this scale with how closely aligned the two are: endpt will have smaller magnitude when they're closely aligned
            other_spatial_steps(k,:) = stepsize*endpt + (stepsize/50)*(-.5+rand(1,3));
        end
        
        other_new_coords = other_old_coords + other_spatial_steps;% + stepsize*(2*(-.5+rand(size(other_old_coords))));
        
        % Now find the nodal branching energy of both the current node, as
        % well as all the othernodes_jth.
        
        new_nodal_branching_energies = zeros((1+length(othernodes_jth)),1);
        
        % First  the current, jth node:
        
        
        
        cf_lengths = zeros(n_cf,1);
        new_unit_vecs = zeros(n_cf,3);
        
        for k=1:n_cf
            cf_lengths(k) = ((other_new_coords(k,1)-newx)^2 + (other_new_coords(k,2)-newy)^2 + (other_new_coords(k,3)-newz)^2)^.5;
            new_unit_vecs(k,1) = (other_new_coords(k,1) - newx)/cf_lengths(k);
            new_unit_vecs(k,2) = (other_new_coords(k,2) - newy)/cf_lengths(k);
            new_unit_vecs(k,3) = (other_new_coords(k,3) - newz)/cf_lengths(k);
        end

        % Now take all the dot products between fibers to see which two are the
        % most aligned.
        N_combinations = nchoosek(n_cf,2);
        %N_combinations = factorial(n_cf)/(factorial(n_cf-2)*2);   % formula for number of fiber combinations at this node
        new_dot_prods = zeros(N_combinations,1);
        combos = combnk(1:n_cf,2);
        for k=1:N_combinations
            new_dot_prods(k) = dot(new_unit_vecs(combos(k,1),:),new_unit_vecs(combos(k,2),:));
        end
        
        [mindot,id]=min(new_dot_prods);
        new_SFs = combos(id,:);
        sf1 = new_SFs(1);   sf2 = new_SFs(2);
        
        new_other_fibers = [1:length(othernodes_jth)]';
        id2 = [sf1;sf2];
        new_other_fibers(id2,:)=[];     % in case there are only 2 fibers at a node, maybe we don't want to put a []
                                % here because that would mess up the
                                % indexing later?
        new_of_stored = new_other_fibers;
        
        term_1 = 1+mindot;
        
        % sf1 and sf2 are the two straightest fibers among the ones at this
        % node. These fibers should be made as straight as possible while all
        % the other fibers need to align towards either sf1 or sf2.
        
        % Dot sf1 with all the others, and then dot sf2 with all the others,
        % and see which one has closer alignment
        
        other_unit_vecs = new_unit_vecs;
        other_unit_vecs([sf1;sf2],:)=[];
        %other_unit_vecs(sf2,:)=[];
        sf1_uv = new_unit_vecs(sf1,:);
        sf2_uv = new_unit_vecs(sf2,:);
        
        dot_prods_sf1 = zeros(size(other_unit_vecs,1),1);
        dot_prods_sf2=dot_prods_sf1;
        
        for k=1:size(other_unit_vecs,1);
            dot_prods_sf1(k) = dot(sf1_uv,other_unit_vecs(k,:));
            dot_prods_sf2(k) = dot(sf2_uv,other_unit_vecs(k,:));
        end
        
        align_sf1 = sum(dot_prods_sf1);
        align_sf2 = sum(dot_prods_sf2);
        
        if align_sf1 > align_sf2
            term_2 = size(other_unit_vecs,1)-align_sf1;
        else
            term_2 = size(other_unit_vecs,1)-align_sf2;
        end
        
        new_nodal_branching_energies(1) = 8*term_1 + term_2;  % New energy of the jth node
        
        % Now we also need to check the branching energies of all the
        % adjacent nodes..........
        
        %other_new_coords = zeros(length(othernodes_jth),3);
        new_other_dot_prods = cell(length(othernodes_jth),1);
        new_all_unit_vecs = cell(length(othernodes_jth),1);
        new_other_SFs = cell(length(othernodes_jth),1);
        new_other_of_stored=new_other_SFs;
        for l=1:length(othernodes_jth)
            current_node = othernodes_jth(l);
            othernodes_lth = othernodes{1,current_node};
            n_cf = length(othernodes_lth);
            new_all_unit_vecs{l}=zeros(n_cf,3);
            cf_lengths = zeros(n_cf,1);
            %other_new_coords(l,:) = nodes(current_node,:);
            newx_1 = other_new_coords(l,1);
            newy_1 = other_new_coords(l,2);
            newz_1 = other_new_coords(l,3);
            
            jth_node_index = find(othernodes_lth==node1);

            for k=1:n_cf
                if k==jth_node_index  % the jth node needs to use the new coordinates, not the ones currently in nodes.
                    cf_lengths(k) = ((newx-newx_1)^2 + (newy-newy_1)^2 + (newz-newz_1)^2)^.5;
                    %unit_vecs(k,:) = (1/cf_lengths(k))*(nodes(othernodes_lth(k,:))-other_new_coords(l,:));
                    new_all_unit_vecs{l}(k,1) = (newx - newx_1)/cf_lengths(k);
                    new_all_unit_vecs{l}(k,2) = (newy - newy_1)/cf_lengths(k);
                    new_all_unit_vecs{l}(k,3) = (newz - newz_1)/cf_lengths(k);
                else
                    %cf_lengths(k) = (sum((nodes(othernodes_lth(k),:)-other_new_coords(l,:)).^2))^.5;
                    cf_lengths(k) = ((nodes(othernodes_lth(k),1)-newx_1)^2 + (nodes(othernodes_lth(k),2)-newy_1)^2 + (nodes(othernodes_lth(k),3)-newz_1)^2)^.5;
                    %unit_vecs(k,:) = (1/cf_lengths(k))*(nodes(othernodes_lth(k,:))-other_new_coords(l,:));
                    new_all_unit_vecs{l}(k,1) = (nodes(othernodes_lth(k),1) - newx_1)/cf_lengths(k);
                    new_all_unit_vecs{l}(k,2) = (nodes(othernodes_lth(k),2) - newy_1)/cf_lengths(k);
                    new_all_unit_vecs{l}(k,3) = (nodes(othernodes_lth(k),3) - newz_1)/cf_lengths(k);
                end
            end
            % Now take all the dot products between fibers to see which two are the
            % most aligned.
            N_combinations = nchoosek(n_cf,2);
            new_other_dot_prods{l}=zeros(N_combinations,1);
            %N_combinations = factorial(n_cf)/(factorial(n_cf-2)*2);   % formula for number of fiber combinations at this node
            combos = combnk(1:n_cf,2);
            for k=1:N_combinations
                new_other_dot_prods{l}(k) = dot(new_all_unit_vecs{l}(combos(k,1),:),new_all_unit_vecs{l}(combos(k,2),:));
            end
            
            [mindot,id]=min(new_other_dot_prods{l});
            new_other_SFs{l} = combos(id,:);
            sf1 = new_other_SFs{l}(1);   sf2 = new_other_SFs{l}(2);
            new_other_fibers = [1:length(othernodes_lth)]';
            id2 = [sf1;sf2];
            new_other_fibers(id2)=[];     % in case there are only 2 fibers at a node, maybe we don't want to put a []
                                % here because that would mess up the
                                % indexing later?
            new_other_of_stored{l} = new_other_fibers;
            
            term_1 = 1+mindot;
            
            % sf1 and sf2 are the two straightest fibers among the ones at this
            % node. These fibers should be made as straight as possible while all
            % the other fibers need to align towards either sf1 or sf2.
            
            % Dot sf1 with all the others, and then dot sf2 with all the others,
            % and see which one has closer alignment
            
            new_other_unit_vecs = new_all_unit_vecs{l};
            new_other_unit_vecs([sf1;sf2],:)=[];
            %other_unit_vecs(sf2,:)=[];
            sf1_uv = new_all_unit_vecs{l}(sf1,:);
            sf2_uv = new_all_unit_vecs{l}(sf2,:);
            
            dot_prods_sf1 = zeros(size(new_other_unit_vecs,1),1);
            dot_prods_sf2 = dot_prods_sf1;
            
            for k=1:size(new_other_unit_vecs,1);
                dot_prods_sf1(k) = dot(sf1_uv,new_other_unit_vecs(k,:));
                dot_prods_sf2(k) = dot(sf2_uv,new_other_unit_vecs(k,:));
            end
            
            align_sf1 = sum(dot_prods_sf1);
            align_sf2 = sum(dot_prods_sf2);
            
            if align_sf1 > align_sf2
                term_2 = size(new_other_unit_vecs,1)-align_sf1;
            else
                term_2 = size(new_other_unit_vecs,1)-align_sf2;
            end
            
            new_nodal_branching_energies(1+l) = 8*term_1 + term_2;
        end
        
        %cnis = [node1;othernodes_jth];
        
        old_nodal_branching_energies = nodal_branching_energy(cnis);
        rp = rand(1);
        
        if sum(new_nodal_branching_energies) < sum(old_nodal_branching_energies)
            %nodes(j,:) = [newx newy newz];  % accept change on current node
            nodes(cnis,:) = [newx newy newz; other_new_coords];
            nodal_branching_energy(cnis) = new_nodal_branching_energies;
            %update stored cells for jth node:
            unit_vecs{j} = new_unit_vecs;
            dot_prods{j} = new_dot_prods;
            SFs{j} = new_SFs;
            of_stored{j} = new_of_stored;
            %update stored cells for othernodes_jth:
            for k=1:length(othernodes_jth)
                unit_vecs{othernodes_jth(k)} = new_all_unit_vecs{k};
                dot_prods{othernodes_jth(k)} = new_other_dot_prods{k};
                SFs{othernodes_jth(k)} = new_other_SFs{k};
                of_stored{othernodes_jth(k)}=new_other_of_stored{k};
            end
            
            N_accepted_BO=N_accepted_BO+1;
        % below is inspired by simulated annealing, helps to jog out of local minima    
        elseif rp<.01 && sum(new_nodal_branching_energies) < 1.5*sum(old_nodal_branching_energies)
            n_baa = n_baa+1;
            nodes(cnis,:) = [newx newy newz; other_new_coords];
            nodal_branching_energy(cnis) = new_nodal_branching_energies;
            %update stored cells for jth node:
            unit_vecs{j} = new_unit_vecs;
            dot_prods{j} = new_dot_prods;
            SFs{j} = new_SFs;
            of_stored{j} = new_of_stored;
            %update stored cells for othernodes_jth:
            for k=1:length(othernodes_jth)
                unit_vecs{othernodes_jth(k)} = new_all_unit_vecs{k};
                dot_prods{othernodes_jth(k)} = new_other_dot_prods{k};
                SFs{othernodes_jth(k)} = new_other_SFs{k};
                of_stored{othernodes_jth(k)}=new_other_of_stored{k};
            end
            
        end
        
    end
    total_BE_intermed(m) = sum(nodal_branching_energy);
    N_accepted_wrt_m(m)=N_accepted_BO;
end

figure
histogram(nodal_branching_energy)
title('Nodal Branching Energy After Branching Optimization')

total_branching_energy_final = sum(nodal_branching_energy)

figure
plot([1:N_branching_optimize],N_accepted_wrt_m)
title('Accepted Changes vs Iteration Number')

figure(40)
plot([1:N_branching_optimize],total_BE_intermed)
title('Branching Energy vs Iteration Number')


end


