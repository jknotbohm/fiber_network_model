% this will randomly delete some fibers as per probability 'prob' of
% existance of a fiber. 
% written by Mainak Sarkar

function [el_set_final] = random_fiber_eliminator(pr, el_set_final)


for k = 1 : size(el_set_final,1)
    if el_set_final(k,1) ~= 0 
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

el_set_final = el_set_final(any(el_set_final,2),:) ;

    

