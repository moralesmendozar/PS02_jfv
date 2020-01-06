function [] = testsmolyak()
           
    %fcn = @(x) 1+sum(min(x.^4,.33),2); %fx = exp(-sum(x.^2,2)); 
    fcn = @(x) sum(x.^3,2); %fx = exp(-sum(x.^2,2)); 
    d = 5;
    mu = 2;
    a = 4*[-1 -1 -1 -1 -1 -1];
    b = 2*[1 1 1 1 1 1];    

    %Get list of points for f to be evaluated at (only has to be done once :)
    [x ienumlist] = smolyakapprox_step1(d,mu,a,b); 
    
    for ind = 1:2000

        fx = fcn(x);

        %Calculate coefficients using calculated values of f
        [ienumlist] = smolyakapprox_step2(fx,ienumlist);

        %Evalutate approximated f at whatever values one wishes
        xpts = 2*rand(32,d)-1; %get random draw from [-1,1]^d
        for dind = 1:d
            xpts(:,dind) = (xpts(:,dind)+1)*(b(dind)-a(dind))/2 + a(dind); %convert to space defined by a,b
        end

        if ind==1
            [fhat,work] = smolyakapprox_step3(ienumlist,xpts);
        else
            [fhat,work] = smolyakapprox_step3(ienumlist,xpts,work);
        end

    end

    ienumlist

    %Compare true value with residual
    plot3(xpts(:,1),xpts(:,2),fhat,'.',xpts(:,1),xpts(:,2),fcn(xpts),'x')    
    %plot(fhat - fcn(xpts))

end


% % Assume a base k number is represented by a k dimensional row vector:
% % [i(1) i(2) ... i(d)] = i(1)*k^d-1 + i(2)*k^d-2 + ... +  i(d-1)*k^1 + i(d)*k^0
% % i should be a dimension k vector.  Use at most n decimal places
% function i10 = convertBaseKtoBaseTen(ik,k)
%     d = size(ik,2);
%     baseconversion = cumprod(repmat(k,[d 1]))/k;
%     baseconversion = baseconversion(d:-1:1);
%     i10=ik*baseconversion;    
% end
% % Assume a base 10 number is represented just like usual.  Want to convert
% % this to a base k number represent by a row vector. 
% % [i(1) i(2) ... i(d)] = i(1)*k^d-1 + i(2)*k^d-2 + ... +  i(d-1)*k^1 + i(d)*k^0
% % Use at most n decimal places.
% function ik = convertBaseTentoBaseK(i10,k,n)
%     ind = 1;
%     prev = i10;
%     
%     baseconversion = cumprod(repmat(k,[n 1]))/k;
%     baseconversion = baseconversion(n:-1:1);
%     % Largest number that can be represented is
%     largenum10 = sum(baseconversion)*k;
%     if any(i10>largenum10) 
%         error('i10 cannot be represented with so few digits')
%     end
%     
%     prev = i10;
%     for ind = n:-1:1
%         ik(:,ind) = rem(prev,k);
%         prev = floor(prev/k);
%     end
%         
% end





