% -------------------------------------------------------------------------
%                           DESCRIPTION
% -------------------------------------------------------------------------
% This function computes a polynomial approximation given a list of values,
% fx, and the structure from smolyak_grid.m, s. If fx is a matrix, then it 
% creates a polynomial approximation treating each column as its own set of
% function values. 
%
% -------------------------------------------------------------------------
%                             SYNTAX
% ------------------------------------------------------------------------- 
%   [s] = smolyak_grid(fx, s)
%
% -------------------------------------------------------------------------
%                               INPUT
% -------------------------------------------------------------------------
%  - s::Structure : results from smolyak_grid.m
%  - fx::Array : 
%
% -------------------------------------------------------------------------
%                               OUTPUT
% -------------------------------------------------------------------------
%  - s::Struct : is a structure containing information used by
%   smolyakapprox_step3.
%
% -------------------------------------------------------------------------
function [s] = smolyak_polynomial(fx,s)
    
    leni = s.leni;

    % Check dims
    if (size(fx,1)~=length(s.redo)) 
        error('fx dims are incorrect')
    end

    % Unpack fx 
    s.f = fx(s.undo,:);

    % Construct the polynomial coefficients
    s.theta = NaN(size(s.f));
    for iind = 1:leni
        for lind = s.Lb(iind):s.Ub(iind)
            tmp = (s.T(s.Lb(iind):s.Ub(iind),lind-s.Lb(iind)+1)./s.cjprod(s.Lb(iind):s.Ub(iind)))';
            s.theta(lind,:) = s.const(iind)/s.clprod(lind)*tmp*s.f(s.Lb(iind):s.Ub(iind),:);
        end 
    end
    
    % Precompute the constant times the coefficient
    s.constTimesTheta = s.theta;
    for iind = 1:leni
        s.constTimesTheta(s.Lb(iind):s.Ub(iind),:) = s.constVec(s.ibar(iind)-s.q+1)*s.constTimesTheta(s.Lb(iind):s.Ub(iind),:);
    end
    
end
