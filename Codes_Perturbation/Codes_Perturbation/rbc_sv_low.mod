// Basic RBC Model with GHH preferences and stochastic volatility
//
// Jesus Fernandez-Villaverde
// Bala Cynwyd, July 9, 2010

//----------------------------------------------------------------
// 0. Housekeeping
//----------------------------------------------------------------

close all
clc

//----------------------------------------------------------------
// 1. Endogenous variables (11=5+4+2)
//----------------------------------------------------------------
 
var 

// Allocation (5)
y c i l k

// Prices (4)
m Rb r w

// Productivity (2)
z ssigmaz;

//----------------------------------------------------------------
// 2. Exogenous variables (2)
//----------------------------------------------------------------
 
varexo 

ez ussigma;

//----------------------------------------------------------------
// 3. Parameters
//----------------------------------------------------------------

parameters 

// Preferences
bbeta ppsi zzeta

// Technology
cte aalpha ddelta rrho ssigmazmean rrhossigmaz

// S.D.'s stochastic processes
eetasigma;

//----------------------------------------------------------------
// 4. Calibration
//----------------------------------------------------------------

// Preferences
bbeta       = 0.99;
zzeta       = 0.5;

// Technology
aalpha      = 1/3;
ddelta      = 0.025;
rrho        = 0.95;  

// Stochastic processes
ssigmazmean = 0.007;
rrhossigmaz = 0.95;
eetasigma   = 0.1;

//----------------------------------------------------------------
// 5. Computation
//----------------------------------------------------------------

// Things we impose
y_ss  = 1;
l_ss  = 1/3;

// Direct consequences
w_ss  = (1-aalpha)*y_ss/l_ss;
ppsi  = w_ss/(l_ss^zzeta);
m_ss  = bbeta;
R_ss  = 1/m_ss;
r_ss  = 1/bbeta-1+ddelta;

// Derivations
k_ss  = aalpha*y_ss/r_ss;
cte   = y_ss/((k_ss^aalpha)*(l_ss^(1-aalpha)));
i_ss  = ddelta*k_ss;
c_ss  = y_ss-i_ss;

//----------------------------------------------------------------
// 6. Model
//----------------------------------------------------------------

model; 
  
  // 1. Pricing kernel
  m = bbeta*((c(-1)-ppsi*(l(-1)^(1+zzeta))/(1+zzeta))/(c-ppsi*(l^(1+zzeta))/(1+zzeta)));

  // 2. Euler equation return on private capital
  1 = m(+1)*(1+r(+1)-ddelta);
  
  // 3. Euler equation return on risk-free bond
  1 = m(+1)*Rb;
  
  // 4. Static condition leisure-consumption
  ppsi*(l^zzeta) = w;

  // 5. Rental rate of capital
  r = aalpha*y/k(-1);

  // 6. Wages
  w = (1-aalpha)*y/l;

  // 7. Law of motion for private capital
  k = (1-ddelta)*k(-1)+i;

  // 8. Production function
  y = cte*((k(-1)^aalpha)*((exp(z)*l)^(1-aalpha)));

  // 9. Resource constraint of the economy 
  c+i = y;

  // 10. Productivity process 
  z = rrho*z(-1)+ssigmazmean*exp(ssigmaz)*ez;

  // 11. Stochastic volatility
  ssigmaz = rrhossigmaz*ssigmaz(-1)+eetasigma*ussigma;

end;

//----------------------------------------------------------------
// 7. Computation
//----------------------------------------------------------------

initval;
  y = y_ss;
  c = c_ss;
  i = i_ss;
  l = l_ss;
  k = k_ss;
  m = m_ss;
  Rb = R_ss;
  r = r_ss;
  w = w_ss;
  z = 0;
  ez = 0;
  ssigmaz = ssigmazmean;
  ussigma = 0;
end;

shocks;
  var ez = 1;
  var ussigma = 1;
end;

steady;

stoch_simul(hp_filter = 1600, irf = 100, order = 2);
