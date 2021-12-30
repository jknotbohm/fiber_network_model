function [inc_nodes,inc_tri] = CreateInclusion(Rmax,Gap)
X=[]; Y=[];

for r=linspace(0,Rmax,ceil(Rmax/Gap))
   if r==0
       xt=0; yt=0;
   else
       arclength = r*2*pi;
       Npts = ceil(arclength/Gap);
       theta = linspace(0,2*pi,Npts);
       theta(end)=[];  %avoid points at both 0 and 2pi
       xt = r*cos(theta);
       yt = r*sin(theta);
   end
   X = [X;xt']; %#ok<*AGROW>
   Y = [Y;yt'];
end
Z=X*0;
inc_nodes = [(1:length(X))',X,Y,Z];
inc_tri = delaunay(X,Y);
inc_tri = [(1:length(inc_tri))',inc_tri];
end

