function [] = testsmolyak()
           
    %fcn = @(x) 1+sum(min(x.^4,.33),2); %fx = exp(-sum(x.^2,2)); 
    fcn = @(x) sum(x.^3,2); %fx = exp(-sum(x.^2,2)); 
    d = 3;
    mu = 2;
    a = -4*ones(1,d);%[-1 -1 -1 -1 -1];
    b = 2*ones(1,d);%*[1 1 1 1 1 1];
    aG = [0.001 0.001 0.001 1 1];
    bG = 2*ones(1,d);

    %Get list of points for f to be evaluated at (only has to be done once :)
    [x ienumlist] = smolyakapprox_step1(d,mu,a,b);
    [numPointsSmolyak d] = size(x);
    
    [GridKT ienumlistKT] = smolyakapprox_step1(d,mu,aG,bG);
    GridKI = GridKT;
    GridB = GridKT;
    GridZ1 = GridKT;
    GridZ2 = GridKT;
    
    ienumlistKI = ienumlistKT;
    ienumlistB = ienumlistKT;
    ienumlistZ1 = ienumlistKT;
    ienumlistZ2 = ienumlistKT;
    
    
    d = fcn(x);  %idea
    options = optimset('Display','off');
    
    d1_11 = zeros(1,numPointsSmolyak);
    d1_12 = zeros(1,numPointsSmolyak);
    d1_21 = zeros(1,numPointsSmolyak);
    d1_22 = zeros(1,numPointsSmolyak);
    d2_11 = zeros(1,numPointsSmolyak);
    d2_12 = zeros(1,numPointsSmolyak);
    d2_21 = zeros(1,numPointsSmolyak);
    d2_22 = zeros(1,numPointsSmolyak);
    d3_11 = zeros(1,numPointsSmolyak);
    d3_12 = zeros(1,numPointsSmolyak);
    d3_21 = zeros(1,numPointsSmolyak);
    d3_22 = zeros(1,numPointsSmolyak);
    d4_11 = zeros(1,numPointsSmolyak);
    d4_12 = zeros(1,numPointsSmolyak);
    d4_21 = zeros(1,numPointsSmolyak);
    d4_22 = zeros(1,numPointsSmolyak);
    
    %find the policy functions in each point of the smolyak grid
    for i=1:numPointsSmolyak
        statesi = GridKT(i,:);
       [ds11 focs11] = fsolve(@(ds)FOCS(ds,statesi,z1,z2,params,dKT11,dKT12,dKT21,dKT22,dKI11,dKI12,dKI21,dKI22),dsinit,options);
%   d1_11(i) = ds(1); 
%   d2_11 = ds(2);
%   d3_11 = ds(3); 
%   d4_11 = ds(4);
    end


    for ind = 1:2000
        
        %Calculate coefficients using calculated values of f
        [ienumlist] = smolyakapprox_step2(d,ienumlist);
        
        [ienumlistKT11] = smolyakapprox_step2(d1_11,ienumlist);
        [ienumlistKI11] = smolyakapprox_step2(d2,ienumlist);
        [ienumlistB11] = smolyakapprox_step2(d3,ienumlist);
%         [ienumlistZ1] = smolyakapprox_step2(d4,ienumlist);
%         [ienumlistZ2] = smolyakapprox_step2(d5,ienumlist);
        
        %Evalutate approximated f at whatever values one wishes
        xpts = 2*rand(32,d)-1; %get random draw from [-1,1]^d
        for dind = 1:d
            xpts(:,dind) = (xpts(:,dind)+1)*(b(dind)-a(dind))/2 + a(dind); %convert to space defined by a,b
        end

        if ind==1
            [fhat,work] = smolyakapprox_step3(ienumlist,xpts);
        else
            [fhat,work] = smolyakapprox_step3(ienumlistKT11,xpts,work);
            
            dKT11 =@(s)smolyakapprox_step3(ienumlistKT11,s,work);
            dKI11 =@(s)smolyakapprox_step3(ienumlistKT11,s,work);
            dKPT11 =@(s)smolyakapprox_step3(ienumlistKT11,s,work);
            dB11 =@(s)smolyakapprox_step3(ienumlistKT11,s,work);
            
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





