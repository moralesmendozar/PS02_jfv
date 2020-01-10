function [y1, extra1, extra2] = ftest(x)
%this function, x^2 is just a test to see if extra1 matters when doing
%   fsolve, solution for x should be +-sqrt(2), depending on inital value,
%           and extra1 should be +-sqrt(2).^3
    y1 = x.^2-2;
    extra1 = x.^2;
    extra2 = x.^3;
    
end