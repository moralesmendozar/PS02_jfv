
function [feval work] = smolyakapprox_step3(s,x,varargin)
%smolyakapprox_step3 Smolyak Polynomial Evaluation
%   [feval] = smolyakapprox_step3(s,x)
%   Given pol which defines an approximation to f, compute the approx
%   value of f at the pts "x".  For performance, it also returns a struct
%   called work that can be used on subsequent calls to improve performance.
% 
%   Grey Gordon 2011.  
%   "But God chose the foolish things of the world to shame the wise; 
%    God chose the weak things of the world to shame the strong." -1Cor12:47

%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.

%
% NOTE: work has been obsoleted (it is now computed in step1 and passed
% along automatically).  
%
% I HAVEN'T IMPLEMENTED THIS FEATURE YET
% If a matrix for fx was given in step2, the coefficients are a matrix defining
% an approximation in each column. If xpts is a matrix, it is then assumed
% all approximations should be evaluated at these points.  If xpts is a
% 3darray with the dim(xpts==3) the same as dim(fx==2), then it is assumed each
% approximation (which is given by a column of .coeffs) should be evaluated by its 
% corresponding value of xpts which are given by (xpts(:,:,columns)).
%

    work = 'The work array will be obsoleted in future versions';

    d = s.d;

    if (d~=size(x,2)), error('smolyakapprox_step3: x is of wrong dimension'),end

    % Convert the x points to z via linear transformation
    z = NaN(size(x));
    for dind = 1:d
        z(:,dind) = 2.d0*(x(:,dind)-s.a(dind))/(s.b(dind)-s.a(dind)) - 1.d0;
    end

    % Evaluate the cheby values beforehand for all orders
    T = ones(size(x,1),s.M_mup1,d);
    %T(:,1,:) = 1.d0
    for dind = 1:d
        T(:,2,dind) = z(:,dind);
    end 
    
    % q<=ibar<=d+mu.  So ibar=d+mu is largest.  
    %Largest any order Cheby can be is the largest any component of m(i) can be.  
    %This occurs when i = 1 except for one component. 
    % => icomp = d+mu-(d-1) = mu+1 => i = mu+1 => m(i) = m(mu+1);
    % NOTE: I precompute this now as s.M_mup1
    twoz = 2*z;
    for oind = 3:s.M_mup1
        for dind = 1:d
            T(:,oind,dind) = twoz(:,dind).*T(:,oind-1,dind) - T(:,oind-2,dind);
        end
    end

    %For each possible value of l, compute the product across i of
    %T(:,li) where li is a component of l    
    Tprod = ones(size(x,1),size(s.l,1));
    for dind = 1:d
        for lind = 1:size(s.l,1)
            if (s.l(lind,dind)>1) 
                Tprod(:,lind) = Tprod(:,lind).*T(:,s.l(lind,dind),dind);
            end
        end
    end

    % This should be much faster
%     tmp_constTimesTheta = NaN(size(s.l,1),1);
%     for iind = 1:leni
%         tmp_constTimesTheta(s.Lb(iind):s.Ub(iind)) = s.constVec(s.ibar(iind)-q+1)*s.theta(s.Lb(iind):s.Ub(iind));
%     end 
    % NOTE: I precompute this constTimesTheta now in step2 
    feval = Tprod*s.constTimesTheta;
    
end
