% Automatic cropper of the 3D fiber network
% written by Mainak Sarkar

function [nodes, fibers, lx, ly, lz] = auto_crop_tool(p)   % p is fraction

index_inp = 1003 ;
load(['C:\Users\msarkar3\Desktop\MATLAB_Path\3D_printing_2D\',...
    '3D_d0p01f30_Network_350x350x5_Seed_',num2str(index_inp),'\'...
    '3D_d0p01f30_Network_350x350x5_Seed_',num2str(index_inp),'.mat'],'nodes','fibers', 'lx', 'ly', 'lz', 'l_fiber') 
% [fibers, nodes] = add_midnode_in_fiber_3D(fibers, nodes) ;
% [fibers, nodes] = add_midnode_in_fiber_3D(fibers, nodes) ;

% lx = (1-p)*max(nodes(:,1)) - (1-p)*min(nodes(:,1)) ;
% ly = (1-p)*max(nodes(:,2)) - (1-p)*min(nodes(:,2)) ;
lx = (1-p)*(lx/2) + (1-p)*(lx/2) ;
ly = (1-p)*(ly/2) + (1-p)*(ly/2) ;

a = find(nodes(:,1) > ((1-p)*max(nodes(:,1)))) ;
b = find(nodes(:,1) < ((1-p)*min(nodes(:,1)))) ;

c = find(nodes(:,2) > ((1-p)*max(nodes(:,2)))) ;
d = find(nodes(:,2) < ((1-p)*min(nodes(:,2)))) ;

% manipulate nodes
nodes = [(1:size(nodes,1))'  nodes(:,1) nodes(:,2) nodes(:,3)] ;


for i = 1 : size(a,1)
    nodes(a(i),:) = NaN ;
end

for i = 1 : size(b,1)
    nodes(b(i),:) = NaN ;
end

for i = 1 : size(c,1)
    nodes(c(i),:) = NaN ;
end

for i = 1 : size(d,1)
    nodes(d(i),:) = NaN ;
end

nodes(~any(~isnan(nodes), 2),:)=[] ;


% operation on fibers:
for j = 1 : size(fibers,1)
    for ii = 1 : length(a)
    if fibers(j,1) == a(ii)
       fibers(j,:) = zeros(1,2) ;
    end
    if fibers(j,2) == a(ii)
       fibers(j,:) = zeros(1,2) ;
    end
    end
end

for j = 1 : size(fibers,1)
    for ii = 1 : length(b)
    if fibers(j,1) == b(ii)
       fibers(j,:) = zeros(1,2) ;
    end
    if fibers(j,2) == b(ii)
       fibers(j,:) = zeros(1,2) ;
    end
    end
end

for j = 1 : size(fibers,1)
    for ii = 1 : length(c)
    if fibers(j,1) == c(ii)
       fibers(j,:) = zeros(1,2) ;
    end
    if fibers(j,2) == c(ii)
       fibers(j,:) = zeros(1,2) ;
    end
    end
end

for j = 1 : size(fibers,1)
    for ii = 1 : length(d)
    if fibers(j,1) == d(ii)
       fibers(j,:) = zeros(1,2) ;
    end
    if fibers(j,2) == d(ii)
       fibers(j,:) = zeros(1,2) ;
    end
    end
end

fibers = fibers(any(fibers,2),:) ;

% renumbering the node numbers and removal of ghost node numbers:
for t = 1 : size(nodes, 1)
  A(t,:) = [t, nodes(t, 1)]  ;
end
    
for uu = 1 : size(fibers,1)
for d = 1 : size(A,1)
if fibers(uu,1) ==  A(d, 2)
   fibers(uu,1) = A(d, 1) ;
end
if fibers(uu,2) ==  A(d, 2)
   fibers(uu,2) = A(d, 1) ;
end
end
end

nodes = [nodes(:,2) nodes(:,3) nodes(:,4)] ;

[fibers, nodes] = check_raw_cropped(fibers, nodes) ;

save(['C:\Users\msarkar3\Desktop\MATLAB_Path\3D_printing_2D\',...
    '3D_d0p01f30_Network_350x350x5_Seed_',num2str(index_inp),'\'...
    '3D_d0p01f30_Network_350x350x5_Seed_cropped_',num2str(index_inp),'.mat'],'nodes','fibers', 'lx', 'ly', 'lz', 'l_fiber') 

