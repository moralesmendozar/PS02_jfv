%UPENN, Econ PhD, 714, prof JFV, student Rodrigo A Morales M
%   Problem set 02, projection method:
function [] = smolyakMain_RAM()
           %This function solves the ps2 model:
           % a DSGE, with 3 endog. states: kT,kI, b
           % and 2 exog states: z1,z2, productivities.
           %    Using the Smolyak Grid and projection
    %% 1. get Parameters:
    params = getParams();
    
    %% 2. SS, get Steady State
    % Solving for [k1_T, k1_I, k_T, K_I, l1, l2] using FOC + LoM
    % FSOLVE W/ DYNARE'S VALUES AS INITIAL GUESS
    x0 =[0.00590718, 0.00191017, 0.00598928,  0.00195713, 0.00393336, 9.32542e-05];
    options = optimoptions('fsolve','Display','iter');
    [x,fval] = fsolve(@(x)ss_foc(x,params),x0,options);

    k1tss = x(1);
    k1iss = x(2);
    kptss = x(3);
    k_Iss = x(4);
    l1ss = x(5);
    l2ss = x(6);

    % Solving for b using FOC w.r.t. b'
    bss = (params.rFirm - params.rBond)/0.02 + 0.2;

    results = table(k1tss,k1iss,kptss,k_Iss,l1ss,l2ss,bss);
    display(results)
    
    %{
                  FSOLVE W/ DYNARE'S VALUES AS INITIAL GUESS
k1_T_ss      k1_I_ss      k_T_ss      k_I_ss        l1_ss        l2_ss       b_ss 
_________    _________    ________    _________    _________    __________    _____
0.0098306    0.0016105    0.010018    0.0016721    0.0045426    0.00010418    -0.05
                  
            DYNARE
value            		 0.115007
k1_T             		 0.00590718
k1_I             		 0.00191017
l1               		 0.00393336
k2_T             		 8.21054e-05
k2_I             		 4.69562e-05
l2               		 9.32542e-05
b                		 -0.05
y1               		 0.0065556
y2               		 0.000195713
s                		 0.00595667
x                		 0.000598928
d                		 0.00280506
k_T              		 0.00598928
k_I              		 0.00195713
z1               		 0
z2               		 0
%}
    
    %% 3. Smolyak
    d = 3;              % number of states
    mu = 2;     % indirect number of Tschebyschow polynomial aproximation
    %a = -4*ones(1,d);   % lower bounds for the states
    aG = [0.001 0.001 -0.8];%0.001];
    a = aG;
    %b = 2*ones(1,d);    % upper bounds for the states
    bG = 0.5*ones(1,d);
    b = bG;

    
    %Get GRID of points for f to b evaluated at (only has 2 be done once :)
    [x, ienumlist] = smolyakapprox_step1(d,mu,a,b);
    numPointsSmolyak = size(x,1); % num of points for smolyak projctn
    nGrid = numPointsSmolyak;
    
    
    % get cartesian product to include the 4 states combination as matrix:
    %    meaning the function will result in four columns (one for each zz)
    exo_states = cartprod(params.z1_grid, params.z2_grid);
    states = cartprod(x,exo_states);
    %states include the z1z2 combinations (each represented by one column)
    nExoStates = size(exo_states,1);   % 4 = # of exog state combinations 
    
    % ------------------
    % Pre-allocation of Policy functions 
    % ------------------
    
    %matrix of functions at points:
    d_k1T = ones(nGrid,4);
    d_k1I = ones(nGrid,4);
    d_kTP = ones(nGrid,4);
    d_b = ones(nGrid,4);
    d_l1 = ones(nGrid,4);
    d_l2 = ones(nGrid,4);
    % group of of functions:
    dKT1 =@(s) k1tss*ones(1,4);
    dKI1 =@(s) k1iss*ones(1,4);
    dKPT =@(s) kptss*ones(1,4);
    dB =@(s) bss*ones(1,4);
    dl1 =@(s) l1ss*ones(1,4);
    dl2 =@(s) l2ss*ones(1,4);
    
    GridKT = x; %rename the grid of points
    ienumlistKT = ienumlist;
    
    
    
    
    
    %get optimal labor1 and 2:
%     z1_vec = states(:,1);
%     d_l1 = labor1(d_k1T,d_k1I,z1_vec,params);
    
    
    
    
    
    % set some parameters for the policy function iteration:
    options = optimset('Display','off');
    tol = 1e-3;
    maxite = 500;
    ite = 0;
    convErr = 10;
    
    % Do Policy Function iteration:
    while( convErr>tol && ite < maxite)
%    for ind = 1:2000
        ite = ite+1;
        
        % Find the policy functions in each point of the smolyak grid
        % solving the FOCs
        if ite ==1
            dsinit = [k1tss,k1iss,kptss, bss];
            for ii=1:numPointsSmolyak
                statesi = GridKT(ii,:);
                for jj = 1:nExoStates
                    zz = exo_states(jj,:);
                    z1 = zz(1);
                    z2 = zz(2);

                    [dssolv, ~, exitfalg] = fsolve(@(ds)FOCS(ds,statesi,z1,z2,params,dKT1,dKI1,dB,dl2),dsinit,options);
                    if exitfalg <1 %if no solution found, keep ss
                        dssolv = dsinit;
                    end
                    %[dssolv focs11] = fsolve(@(ds)FOCS(ds,statesi,params,dKT1,dKI1),dsinit,options);
                    [~, l1d, l2d]  = FOCS(dssolv,statesi,z1,z2,params,dKT1,dKI1,dB,dl2);
                    d_k1T(ii,jj)    = dssolv(1);
                    d_k1I(ii,jj)    = dssolv(2);
                    d_kTP(ii,jj)    = dssolv(3);
                    d_b(ii,jj)      = dssolv(4);
                    d_l1(ii,jj)     = l1d;
                    d_l2(ii,jj)     = l2d;

                end %for states zz ( solving policy functions at smolyak
            end  %for ii policy functions
        else
            for ii=1:numPointsSmolyak
                statesi = GridKT(ii,:);
                for jj = 1:nExoStates
                    zz = exo_states(jj,:);
                    z1 = zz(1);
                    z2 = zz(2);
                    dsinit = [dKT1previous(ii,jj),dKI1previous(ii,jj),dKPTprevious(ii,jj), dBprevious(ii,jj)];
                    [dssolv, ~, exitfalg] = fsolve(@(ds)FOCS(ds,statesi,z1,z2,params,dKT1,dKI1,dB,dl2),dsinit,options);
                    if exitfalg <1   %if no solution found, keep previous solution
                        dssolv = dsinit;
                    end
                    %[dssolv focs11] = fsolve(@(ds)FOCS(ds,statesi,params,dKT1,dKI1),dsinit,options);
                    [~, l1d, l2d]  = FOCS(dssolv,statesi,z1,z2,params,dKT1,dKI1,dB,dl2);
                    d_k1T(ii,jj)    = dssolv(1);
                    d_k1I(ii,jj)    = dssolv(2);
                    d_kTP(ii,jj)    = dssolv(3);
                    d_b(ii,jj)      = dssolv(4);
                    d_l1(ii,jj)     = l1d;
                    d_l2(ii,jj)     = l2d;

                end %for states zz ( solving policy functions at smolyak
            end  %for ii policy functions
        end
        
        %Calculate coefficients using calculated values of decision policis       
        [ienumlistKT1] = smolyakapprox_step2(d_k1T,ienumlist);
        [ienumlistKI1] = smolyakapprox_step2(d_k1I,ienumlist);
        [ienumlistKTp] = smolyakapprox_step2(d_kTP,ienumlist);
        [ienumlistBp] = smolyakapprox_step2(d_b,ienumlist);
        [ienumlistL1] = smolyakapprox_step2(d_l2,ienumlist);
        [ienumlistL2] = smolyakapprox_step2(d_l2,ienumlist);
        
        %Evalutate approximated decision at random points
        num_xpts = 500;
        xpts = 2*rand(num_xpts,d)-1; %get random draw from [-1,1]^d
        %convert to space defined by a,b
        for dind = 1:d
            xpts(:,dind) = (xpts(:,dind)+1)*(b(dind)-a(dind))/2 + a(dind); 
        end
        
        dKT1previous = dKT1;
        dKI1previous = dKI1;
        dKPTprevious = dKPT;
        dBprevious = dB;
        dl1previous = dl1;
        dl2previous = dl2;
        dKT1 = @(s)smolyakapprox_step3(ienumlistKT1,s);
        dKI1 = @(s)smolyakapprox_step3(ienumlistKI1,s);
        dKPT = @(s)smolyakapprox_step3(ienumlistKTp,s);
        dB = @(s)smolyakapprox_step3(ienumlistBp,s);
        dl1 =@(s)smolyakapprox_step3(ienumlistL1,s);
        dl2 =@(s)smolyakapprox_step3(ienumlistL2,s);
        
        %%% Compute the convergence error:
        convErr = norm(dKT1(xpts)-dKT1previous(xpts),2);
        
        if( mod(ite,10)== 0 || ite ==1)
            display(['iteration = ', num2str(ite), '. And error = ', num2str(convErr)]);
        end %if ite for display
    end %policy iteration while
    ite
    %ienumlist

    %Compare true value with residual
    matrix_dKT1 = dKT1(xpts); %result is matrix of size:(n_xpts x 4)
    dKT1_11 = matrix_dKT1(:,1);  %get policy function for lowest states
    
    %plot for only two dimensions of xpts (xpts has 3 dimensions!!):
    plot3(xpts(:,1),xpts(:,2),dKT1_11,'.')%,xpts(:,1),xpts(:,2),dKT1_11,'x')    
    %plot(fhat - fcn(xpts))

end