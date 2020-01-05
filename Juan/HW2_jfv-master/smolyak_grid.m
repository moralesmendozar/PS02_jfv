% -------------------------------------------------------------------------
%                           DESCRIPTION
% -------------------------------------------------------------------------
% This function constructs the Smolyak grid points for a hypercube in R^d
% defined by bounds a,b with level of approximation mu
%
% -------------------------------------------------------------------------
%                             SYNTAX
% ------------------------------------------------------------------------- 
%   [x s] = smolyak_grid(d,mu,a,b)
%
% -------------------------------------------------------------------------
%                               INPUT
% -------------------------------------------------------------------------
%  - d::Integer : # of state variables
%  - mu::Integer : level of the approximation 
%  - a::Integer :
%  - b::Integer : 
%
% -------------------------------------------------------------------------
%                               OUTPUT
% -------------------------------------------------------------------------
%  - x::Array : collocation points that your function should be evaluated.
%  - s::Struct : is a structure containing information used by
%   smolyakapprox_step2 and smolyakapprox_step3.
%
% -------------------------------------------------------------------------
function [x, s] = smolyak_grid(d,mu,a,b)

    % 0. Housekeeping
    q = max(d,mu+1);        
    
    s.q = q;
    s.d = d;
    s.mu = mu;
    s.a = a(1:d);
    s.b = b(1:d);
    s.M_mup1 = smolyakM(mu+1); %This constant is used in step3
    
    if (length(a)>d), warning('length of a is longer than dim d'), end
    if (length(b)>d), warning('length of b is longer than dim d'), end
    
    % 1. Construct necessary coefficients using fx values.
        % First, enumerate all necessary theta
        % - For all i satisfying q <= ibar <= d + mu, 
        % - need { l | l(1) is in 1...m(i1), l(2) is in 1...m(i1)}
    
    % Enumerate all such i, then all such m(i), then all such l
    
    % Determine all i that fit the criterion q<=|i|<=d+mu
    tmp=[];
    for ibar = q:d+mu
        enum = smolyakEnumerate(d,ibar-d);        
        enum = enum+1;   
        tmp = [enum;tmp];
    end

    tmp = unique(tmp,'rows');
    leni = size(tmp,1);

    s.leni = leni;
    s.i = tmp;
    s.ibar = sum(s.i,2);
    s.k = smolyakM(s.i);
    
  
    s.Lb = NaN(s.leni,1);
    s.Ub = NaN(s.leni,1);    
    tmpInt = 0;
    for iind = 1:leni        
        s.Lb(iind) = tmpInt + 1;
        tmpInt = tmpInt + prod(s.k(iind,:));
        s.Ub(iind) = tmpInt;
    end
    

    s.z = NaN(s.Ub(leni),d);
    s.j = NaN(s.Ub(leni),d);
    
    for iind = 1:leni
        ztmp = [];
        jtmp = [];
        for dind = 1:d
            ztmp = cartprod(ztmp,smolyakG(s.k(iind,dind)));
            jtmp = cartprod(jtmp,(1:s.k(iind,dind))');
        end
        s.z(s.Lb(iind):s.Ub(iind),:) = ztmp;
        s.j(s.Lb(iind):s.Ub(iind),:) = jtmp;
    end
    
   
    for dind = 1:d
        s.x(:,dind) = (s.z(:,dind)+1.d0)*(b(dind)-a(dind))/2.d0 + a(dind);
    end    
    
    s.f = NaN(size(s.z));
    s.l = s.j; % l and j have the same form but are used differently
    
      
    s.T = NaN(s.Ub(leni),max(s.Ub-s.Lb)+1);
    for iind = 1:leni
        for lind = s.Lb(iind):s.Ub(iind)
            for jind = s.Lb(iind):s.Ub(iind)
                tmpProd = 1.d0;
                for dind = 1:d
                    tmpProd = tmpProd*ChebEvalOneDim(s.l(lind,dind)-1,s.z(jind,dind));
                end
                s.T(jind,lind-s.Lb(iind)+1) = tmpProd;
            end
        end       
    end

    s.clprod = NaN(size(s.z,1),1);
    for iind = 1:leni
        for lind = s.Lb(iind):s.Ub(iind)
            s.clprod(lind) = 1.d0;
            for dind = 1:d
                % Dimensions with only one dimensions are "dropped" according to Krueger and Kubler
                if (s.k(iind,dind)>1) 
                    % The product is 2^(l==1 or l==k)
                    if ((s.l(lind,dind)==1) || (s.l(lind,dind)==s.k(iind,dind))) 
                        s.clprod(lind) = s.clprod(lind)*2.d0;
                    end 
                end 
            end 
        end 
    end 

    % cjprod is used differently but has the same values as clprod 
    s.cjprod = s.clprod;

    % for each i, there is a unique "const"
    s.const = NaN(leni,1);
    for iind = 1:leni
        tmpProd = 1.d0;
        tmpSum = 0.d0;
        for dind = 1:d
            if (s.k(iind,dind)>1) 
                tmpSum = tmpSum + 1.d0;
                tmpProd = tmpProd*(s.k(iind,dind)-1);
            end
        end
        s.const(iind) = (2.d0^tmpSum)/tmpProd;
    end
    
    % Store the "choose" constant values
    s.constVec = NaN(d+mu+1-q,1);
    for ibar = q:d+mu
        s.constVec(ibar-q+1) = (-1.d0)^(d+mu-ibar)*choose(d-1,d+mu-ibar);
    end

    % Get the smolyak grid points (b/c of their nested nature, only return unique values)
    [x s.redo s.undo] = unique(s.x,'rows');
    
end