// Simple New Keynesian Model
// See nke_description.tex
//
// Jesus Fernandez-Villaverde
// Bala Cynwyd, June 26, 2009

//----------------------------------------------------------------
// 0. Housekeeping
//----------------------------------------------------------------

close all

//----------------------------------------------------------------
// 1. Endogenous variables (19=1+5+7+1+2+3)
//----------------------------------------------------------------
 
var 

// Utility (1)
lagrangian

// Allocation (5)
y c x l k

// Prices (7)
m r w mc ppi ppistar Rn

// Inefficiencies (1)
v

// Auxiliary functions (2)
g1 g2

// Stochastic processes (3)
d phi z;

//----------------------------------------------------------------
// 2. Exogenous variables (4)
//----------------------------------------------------------------
 
varexo 

ed ephi ez em;

//----------------------------------------------------------------
// 3. Parameters
//----------------------------------------------------------------

parameters 

// Preferences
bbeta ppsi zzeta

// Technology
aalpha A ddelta eepsilon

// Nominal rigidities
ttheta

// Taylor rule
bigppi ggamma

// Autoregressive stochastic processes
rrhod rrhophi rrhoz

// S.D.'s stochastic processes
ssigmad ssigmaphi ssigmaz ssigmam;

//----------------------------------------------------------------
// 4. Calibration
//----------------------------------------------------------------

// Preferences
bbeta     = 0.99;
zzeta     = 0.5;

// Technology
aalpha    = 0.33;
ddelta    = 0.025;
eepsilon  = 10;

// Nominal rigidities
ttheta    = 0.75;

// Taylor rule
bigppi    = 1.005;
ggamma    = 1.5;

// Autoregressive stochastic processes
rrhod     = 0.95;
rrhophi   = 0.95;
rrhoz     = 0.95;

// S.D.'s stochastic processes
ssigmad   = 0.025;
ssigmaphi = 0.01;
ssigmaz   = 0.007;
ssigmam   = 0.001;

// Two scale parameters: A and ppsi
y_ss  = 1;
l_ss  = 1/3;

ppistar_ss = ((1-ttheta*(bigppi^(eepsilon-1)))/(1-ttheta))^(1/(1-eepsilon));
mc_ss      = ((eepsilon-1)/eepsilon)*ppistar_ss*(1-bbeta*ttheta*(bigppi^eepsilon))/(1-bbeta*ttheta*(bigppi^(eepsilon-1)));
v_ss       = (1-ttheta)/(1-ttheta*(bigppi^eepsilon))*(ppistar_ss^(-eepsilon));
r_ss       = 1/bbeta-1+ddelta;
w_ss       = (1-aalpha)*v_ss*mc_ss/l_ss;
k_ss       = (aalpha/(1-aalpha))*w_ss/r_ss*l_ss;
x_ss       = ddelta*k_ss;
A          = (v_ss)/((k_ss^aalpha)*(l_ss^(1-aalpha)));
c_ss       = y_ss-ddelta*k_ss;
lag_ss     = 1/c_ss;
ppsi       = w_ss/(c_ss*(l_ss^zzeta));

//----------------------------------------------------------------
// 5. Model
//----------------------------------------------------------------

model; 
   
  // 1. Lagrangian
  lagrangian = exp(d)/c;

  // 2. Pricing kernel
  m = bbeta*lagrangian/lagrangian(-1);

  // 3. Euler equation return on private capital
  1 = m(+1)*(r(+1)+1-ddelta);
  
  // 4. Euler equation return on bonds
  1 = m(+1)*Rn/ppi(+1);

  // 5. Leisure choice
  exp(phi)*ppsi*(l^zzeta) = w/c;  

  // 6. Recursive equation for prices 1
  eepsilon*g1 = (eepsilon-1)*g2;

  // 7. Recursive equation for prices 2 
  g1 = lagrangian*mc*y+bbeta*ttheta*(ppi(+1)^eepsilon)*g1(+1);

  // 8. Recursive equation for prices 3
  g2 = lagrangian*ppistar*y+bbeta*ttheta*(ppi(+1)^(eepsilon-1))*(ppistar/ppistar(+1))*g2(+1);

  // 9. FOC of firms with respect to capital and labor
  k(-1)/l = aalpha/(1-aalpha)*w/r;

  // 10. Marginal cost
  mc = ((1/(1-aalpha))^(1-aalpha))*((1/aalpha)^aalpha)*(w^(1-aalpha))*(r^aalpha)/(A*exp(z));
      
  // 11. Evolution of prices
  1 = ttheta*(ppi)^(eepsilon-1)+(1-ttheta)*ppistar^(1-eepsilon);
  
  // 12. Taylor rule
  Rn = (1/bbeta)*(bigppi^(1-ggamma))*(ppi^ggamma)*exp(ssigmam*em);

  // 13. Production function
  y = (A/v)*(k(-1)^aalpha)*((exp(z)*l)^(1-aalpha));

  // 14. Resource constraint of the economy 
  c+x = y;

  // 15. Law of motion for private capital
  k = (1-ddelta)*k(-1)+x;

  // 16. Evolution of price dispersion inefficiency
  v = ttheta*(ppi^eepsilon)*v(-1)+(1-ttheta)*(ppistar^(-eepsilon));

  // 17. Intertemporal shock
  d = rrhod*d(-1)+ssigmad*ed;

  // 18. Intratemporal shock
  phi = rrhophi*phi(-1)+ssigmaphi*ephi;

  // 19. Productivity process 
  z = rrhoz*z(-1)+ssigmaz*ez;
  
end;

//----------------------------------------------------------------
// 6. Computation
//----------------------------------------------------------------

initval;
  lagrangian = lag_ss;
  m    = bbeta;
  Rn   = bigppi/bbeta;
  r    = r_ss;
  ppi  = bigppi;
  ppistar = ppistar_ss;
  v    = v_ss;
  mc   = mc_ss;
  w    = w_ss;
  k    = k_ss;
  l    = l_ss;
  x    = ddelta*k_ss;
  y    = y_ss;
  c    = c_ss;
  d    = 0;
  phi  = 0;
  z    = 0;
  ed   = 0;
  ephi = 0;
  ez   = 0;
  em   = 0;
end;

shocks;
    // var VARIABLE_NAME = EXPRESSION;
    // Specifies the variance of a variable.
  var ed   = 1;   
  var ephi = 1;
  var ez   = 1;
  var em   = 1;
end;

steady;

stoch_simul(hp_filter = 1600, irf = 20, order = 1) y c x l k w r ppi Rn;
