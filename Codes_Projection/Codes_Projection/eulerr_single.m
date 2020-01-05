function [value,kp,c,l,y,rb,rk,rkCond,euler_error] = ... 
        eulerr_single(alpha,beta,delta,gamma,nu,rho,...
        Z,PI,k_min,k_max,node_num,shock_num,M,z_index,k)

    z = Z(z_index);

    k_scaled = 2*(k-k_min)/(k_max - k_min) -1;

    Tk = zeros(node_num,1);
    for i = 1:node_num % ith polynomial
        Tk(i) = cos(real(i-1)*acos(k_scaled));
    end

    rho1 = rho(1:M,1);     % coefficient for value function
    rho2 = rho(M+1:2*M,1); % coefficient for labor    
    rho1_section = rho1(((z_index-1)*node_num+1):z_index*node_num);
    rho2_section = rho2(((z_index-1)*node_num+1):z_index*node_num);

    value = dot(rho1_section,Tk);   % Value function at each collocation points    
    l = dot(rho2_section,Tk);       % Labor at each collocation points

    y = exp(z)*k^alpha*l^(1-alpha);
    c = nu/(1-nu)*(1-alpha)*y/l*(1-l);            
    kp = y+(1-delta)*k-c;

    g_k_scaled = 2*(kp-k_min)/(k_max - k_min) -1;    
    T_g_k = zeros(node_num,1);

    for i = 1:node_num % ith polynomial
        T_g_k(i) = cos(real(i-1)*acos(g_k_scaled));
    end

    % calculate residual  
    vp = zeros(shock_num,1);
    for zp_index = 1:shock_num
        rho1_section = rho1(((zp_index-1)*node_num+1):zp_index*node_num);
        vp(zp_index) = dot(rho1_section,T_g_k);
    end

    mp   = zeros(shock_num,1);
    fkp  = zeros(shock_num,1);
    temp = zeros(shock_num,1);

    for zp_index = 1:shock_num

        rho2_section = rho2(((zp_index-1)*node_num+1):zp_index*node_num);
        lp = dot(rho2_section,T_g_k);

        yp = exp(Z(zp_index))*kp^alpha*lp^(1-alpha);
        cp = nu/(1-nu)*(1-alpha)*yp/lp*(1-lp);

        product1 = (cp^nu*(1-lp)^(1-nu))/(c^nu*(1-l)^(1-nu));
        product2 = c/cp;
        mp(zp_index) = beta*product1^((1-gamma))*product2;

        Ucp = nu*((cp^nu*(1-lp)^(1-nu))^(1-gamma))/cp;
        fkp(zp_index) = alpha*exp(Z(zp_index))*kp^(alpha-1)*lp^(1-alpha);
        temp(zp_index) = Ucp*(fkp(zp_index)+1-delta);
    
    end    

    euler_rhs = beta*dot(PI(z_index,:),temp);

    A = euler_rhs/(nu*(1-l)^((1-nu)*(1-gamma)));
    euler_error = 1- A^(1/(nu*(1-gamma)-1))/c;

    euler_error = log10(abs(euler_error));

    rb = (dot(PI(z_index,:),mp))^(-1)-1;    % bond return between t and t+1
    rk = alpha*y/k-delta;                   % risky asset return between t-1 and t
    rkCond = dot(PI(z_index,:),fkp);        % conditional rk next period

    if( euler_error < -17 )
        euler_error = -17;
    end