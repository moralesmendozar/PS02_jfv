% Given the dimensionality of the state, n, enumerates all possible ways to 
% add one mu times to n (where adding one can occur in any of d positions).
% q = n + mu
function [enum] = smolyakEnumerate(n,mu)
    
    if (mu==0)
        enum = zeros([1 n]);
        return
    end
    
    if (mu>=1) 
       
        enum_mum1 = smolyakEnumerate(n,mu-1);
        m = size(enum_mum1,1);
        
        % Given all previous enumerations, I can add one in n places.
        for i = 1:n
            enum(1+(i-1)*m:i*m,:) = enum_mum1;            
            enum(1+(i-1)*m:i*m,i) = enum(1+(i-1)*m:i*m,i) + 1;
        end
        
        return
    end
    
end