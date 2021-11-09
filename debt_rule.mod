var y k c n w r rk i g mc_r pi a q t b;  //Endogenous variables & AR processes
varexo epsilon_g;              //Exogenous shocks 

parameters alpha beta delta varphi phi phi_pi sigmaC thetaN thetaT Gamma rhog rhoa sigma theta eta lambda phi_b phi_g
gy_ss r_ss cy_ss iy_ss;  // Parameters

//--- Parameterization ---//
alpha=1/3;
beta = .99;
delta = .025;
varphi = 0.2;
rhog=.9;
rhoa=0.9;
sigma=0.01;
theta=0.75;
phi=6;
phi_pi=1.5;
eta=1;
lambda=0.5;
phi_b=0.33;
phi_g=0.1;

//--- Steady state expressions ---//
gy_ss=.20;
r_ss=(1/beta)-1;
mc_ss=(phi-1)/phi;
cy_ss= 1-gy_ss-((delta*alpha)/((r_ss+delta)/mc_ss));
iy_ss=1-cy_ss-gy_ss;
Gamma=1/((1/mc_ss)*varphi*cy_ss+(1-alpha)*(1-lambda*(1+varphi))); 
sigmaC=(1-lambda)*Gamma*((1/mc_ss)*varphi*cy_ss+(1-alpha));
thetaN=lambda*Gamma*(1-alpha)*(1+varphi)*varphi; 
thetaT=lambda*Gamma*(1/mc_ss)*varphi; 

//--- Equations of the model (in log-linear approximation) ---//
model(linear);

    %% Household
	% Euler Ricardian
	c = c(+1) -sigmaC*(r-pi(+1))-thetaN*(n(+1)-n)+thetaT*(t(+1)-t);
	% Labour supply
	w = c + varphi*n;
    % Tobin's Q
    q=beta*q(+1)+(1-beta*(1-delta))*rk(+1)-(r-pi(+1)); 
    % Investment
    i-k(-1)=eta*q;
	% Capital law of motion
	k = (1-delta)*k(-1) + delta*i;

	% Intermediary firms
	% Production function
	y = a+ alpha*k(-1) + (1-alpha)*n;
	% Real marginal cost
	mc_r=-a +alpha*rk+(1-alpha)*w;
	% Cost minimization
	w + n = rk + k(-1);
	% Price dynamics NKPC
	pi = beta*pi(+1) + (((1-theta)*(1-theta*beta))/theta)*mc_r;

	%%Government
	%Budget constraint
	b =(1 + r_ss)*(b(-1) + g - t);
	%Fiscal policy rule
    t = phi_b*b(-1) + phi_g*g;

	% Market Clearing
	y = cy_ss*c + iy_ss*i + g;

	% monetary policy
	r = phi_pi*pi;

	% Exogenous shocks
	g = rhog*g(-1) + epsilon_g;
    a = a(-1);

end;
check;
steady;

//--- Shocks ---//
shocks;
var epsilon_g; stderr sigma;
end;

//--- Simulation ---//
stoch_simul(order=1, irf=20, periods=1000);