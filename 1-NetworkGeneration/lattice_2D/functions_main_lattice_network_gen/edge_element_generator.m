% finding neighbors for edge nodes
% written by Mainak Sarkar


function [el_set_1, el_set_2, el_set_3, el_set_4, el_set_5, el_set_6, el_set_7, el_set_8] = edge_element_generator(Nx, Ny, N1)



for m = 1  % first row
     for n = 1:Nx
         if n == 1 
           el_set_1 = [N1(1,1), N1(2,1) ; N1(1,1), N1((Nx+1),1)] ;
         elseif n == Nx
                      if rem(Nx,2) == 0
           el_set_2 = [N1(Nx,1), N1((Nx-1),1); N1(Nx,1), N1(2*Nx,1)] ;
                      end
                      if rem(Nx,2) ~= 0
           el_set_2 = [N1(Nx,1), N1((Nx-1),1); N1(Nx,1), N1(2*Nx,1) ; N1(Nx,1), N1(2*Nx-1,1)] ;
                      end
         end
         
     end
end

i = 0 ; 
if rem(Nx,2) ~= 0 || rem(Nx,2) == 0
for in = 2:(Nx-1)
    i = i + 1 ;
    el_set_5{i} = [N1(in,1), N1(in-1,1); N1(in,1), N1(in+1,1); N1(in,1), N1(Nx+in-1,1) ; N1(in,1), N1(Nx+in,1)] ;
end   
end

el_set_5 = cell2mat(el_set_5') ;

%{
if rem(n,2)==0
    output=even;
else
    output=odd;
end
%}

for m = Ny % last row at the top
if rem(Ny,2)~=0 % for odd last row
    for n = 1:Nx
         if n == 1 
           el_set_3 = [N1(Nx*(Ny-1)+1,1), N1(Nx*(Ny-1)+2,1) ; N1(Nx*(Ny-1)+1,1), N1(Nx*(Ny-2)+1,1)] ; % top left
         elseif n == Nx
           el_set_4 = [N1(Nx*Ny,1), N1((Nx*Ny-1),1); N1(Nx*Ny,1), N1(Nx*(Ny-1),1)] ; % top right
         end
    end
elseif rem(Ny,2) == 0
    for n = 1:Nx
         if n == 1 
           el_set_3 = [N1(Nx*(Ny-1)+1,1), N1(Nx*(Ny-1)+2,1) ; N1(Nx*(Ny-1)+1,1), N1(Nx*(Ny-2)+1,1); N1(Nx*(Ny-1)+1,1), N1(Nx*(Ny-2)+2,1)] ; % top left
         elseif n == Nx
           el_set_4 = [N1(Nx*Ny,1), N1((Nx*Ny-1),1); N1(Nx*Ny,1), N1(Nx*(Ny-1),1) ; N1(Nx*Ny,1), N1(Nx*(Ny-1)-1,1)] ; % top right
         end
    end
end
end

i = 0 ;
for in = 2 : (Nx-1)
    i = i + 1 ;
    el_set_6{i} = [N1(Nx*(Ny-1)+in,1), N1(Nx*(Ny-1)+in-1,1); 
        N1(Nx*(Ny-1)+in,1), N1(Nx*(Ny-1)+in+1,1); 
        N1(Nx*(Ny-1)+in,1), N1(Nx*(Ny-2)+in,1) ; 
        N1(Nx*(Ny-1)+in,1), N1(Nx*(Ny-2)+in+1,1)] ;
end

el_set_6 = cell2mat(el_set_6') ;

% left and right columns
i = 0 ; 
for n = 1
    if rem(Ny,2)~=0
    for m = 3:2:Ny-2
        i = i + 1 ;
        el_set_7{i} = [N1(Nx*(m-1)+1,1), N1(Nx*(m-1)+2,1); N1(Nx*(m-1)+1,1), N1(Nx*(m)+1,1) ; N1(Nx*(m-1)+1,1), N1(Nx*(m-2)+1,1)] ;
    end
    elseif rem(Ny,2) == 0 
    for m = 3:2:(Ny-1)
        i = i + 1 ;
        el_set_7{i} = [N1(Nx*(m-1)+1,1), N1(Nx*(m-1)+2,1); N1(Nx*(m-1)+1,1), N1(Nx*(m)+1,1) ; N1(Nx*(m-1)+1,1), N1(Nx*(m-2)+1,1)] ;
    end
    end
end

el_set_7 = cell2mat(el_set_7') ;

   
i = 0 ; 
for n = Nx
    if rem(Nx,2) == 0 
    if rem(Ny,2)~=0
    for m = 3:2:Ny-2
        i = i + 1 ;
        el_set_8{i} = [N1(Nx*m,1), N1(Nx*m-1,1); N1(Nx*m,1), N1(Nx*(m+1),1) ; N1(Nx*m,1), N1(Nx*(m-1),1)] ;
    end
    elseif rem(Ny,2) == 0 
    for m = 3:2:(Ny-1)
        i = i + 1 ;
        el_set_8{i} = [N1(Nx*m,1), N1(Nx*m-1,1); N1(Nx*m,1), N1(Nx*(m+1),1) ; N1(Nx*m,1), N1(Nx*(m-1),1)] ;
    end
    end
    end
end

i = 0 ; 
for n = Nx
    if rem(Nx,2) ~= 0 
    if rem(Ny,2)==0
    for m = 2:2:(Ny-2)
        i = i + 1 ;
el_set_8{i} = [N1(Nx*m,1), N1(Nx*m-1,1); N1(Nx*m,1), N1(Nx*(m+1),1) ; N1(Nx*m,1), N1(Nx*(m-1),1)] ;
    end
    elseif rem(Ny,2) ~= 0 
    for m = 2:2:(Ny-1)
        i = i + 1 ;
        el_set_8{i} = [N1(Nx*m,1), N1(Nx*m-1,1); N1(Nx*m,1), N1(Nx*(m+1),1) ; N1(Nx*m,1), N1(Nx*(m-1),1)] ;
    end
    end
    end
end


el_set_8 = cell2mat(el_set_8') ;

