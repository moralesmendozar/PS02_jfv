% Given an integer i, deliver the function m(i) = 2^(i-1) + 1
function order = smolyakM(i)        
    
    order = NaN(size(i));
    
    if sum(i<1)>0
        error('i must be >= 1')
    else
        order(i==1) = 1;
        order(i>1) =  2.^(i(i>1)-1) + 1;
    end
    
end