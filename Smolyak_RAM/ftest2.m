function [y1] = ftest2(x1,x2)
%this function, x^2 is just a test to see if extra1 matters when doing
%   fsolve, solution for x should be +-sqrt(2), depending on inital value,
%           and extra1 should be +-sqrt(2).^3
    y1 = x1.^2 +x2.^2;
    
end