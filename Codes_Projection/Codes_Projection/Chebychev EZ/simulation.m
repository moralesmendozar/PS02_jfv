function [vSeriesDrop,kSeriesDrop,cSeriesDrop,lSeriesDrop,ySeriesDrop,...
    eeSeriesDrop,rbSeriesDrop,rkSeriesDrop,rkCondSeriesDrop] = ...
    simulation(alpha,beta,delta,gamma,theta,nu,kss,rho,Z,PI,k_min,k_max,node_num,shock_num,M,T,dropT)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation of the economy for T periods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vSeries=zeros(T,1); 
kSeries=zeros(T,1);
cSeries=zeros(T,1);
lSeries=zeros(T,1);
ySeries=zeros(T,1);
eeSeries=zeros(T,1);
rbSeries=zeros(T,1);
rkSeries=zeros(T,1);
rkCondSeries=zeros(T,1);

% Generate Markov Chain
rand('state',0);
mu=zeros(1,shock_num);      % initial distribution
mu(1,round(shock_num/2))=1; % Always start from Z=0
z_index=zeros(1,T+1);   % Stores indexes of the realization of "z" over
                    % the simulation
z_index(1)=rando(mu); % generate first ind value (time 0, not time 1)

for i=1:T,
    z_index(i+1) = rando(PI(z_index(i),:));
end
kSeries(1) = kss;
for t_index = 1:T;
    [vSeries(t_index),kSeries(t_index+1),cSeries(t_index),lSeries(t_index),...
        ySeries(t_index),rbSeries(t_index),...
        rkSeries(t_index),rkCondSeries(t_index),eeSeries(t_index)] = ...
        eulerr_single(alpha,beta,delta,gamma,theta, ...
        nu,rho,Z,PI,k_min,k_max,node_num,shock_num,M,z_index(t_index),kSeries(t_index));
end

vSeriesDrop = vSeries(dropT+1:T); 
kSeriesDrop = kSeries(dropT+1:T);    
cSeriesDrop = cSeries(dropT+1:T);
lSeriesDrop = lSeries(dropT+1:T);
ySeriesDrop = ySeries(dropT+1:T);
eeSeriesDrop = eeSeries(dropT+1:T);
rbSeriesDrop = rbSeries(dropT+1:T);
rkSeriesDrop = rkSeries(dropT+1:T);
rkCondSeriesDrop = rkCondSeries(dropT+1:T);

