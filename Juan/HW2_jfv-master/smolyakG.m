% Given an integer m, construct the set of m extrema of the Chebyshev 
% polynomials x_j = -cos(pi*(j-1)/(n-1)), j = 1,...,m
function [set j]= smolyakG(m)
    
    if (m<1) 
        
        error('smolyakG: m must be >=1')
    
    elseif (m==1)
    
        set = 0;
        j = (1:m)';
    
    else
        
        j = (1:m)';
        set = -cos(pi*(j-1)/(m-1));
        set(abs(set)<1e-12) = 0;
        set(abs(set-1.e0)<1e-12) = 1;
        set(abs(set+1.e0)<1e-12) = -1;
    
    end
   
end