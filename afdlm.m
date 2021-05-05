function [y_pred, theta, lambda, e_t, dedl, step] = afdlm(X, y, theta0, S0, lambda0, lambda_max0, lambda_min0, beta0,a,b1,b2)
%ADAPTIVE FORGETTING DYNAMIC LINEAR MODEL 
% [Y_PRED, THETA] = AFDLM(X, Y, L0, LMAX, LMIN, B, A, B1, B2) 
% implements the Adaptive Forgetting Dynamic Linear Model algorithm using the stochastic
% gradient descent algorithm ADAM to incrementally tune the forgetting factor (\lambda)
% 

% Inputs:
%   -- DATA --
%	(X): T-by-D data matrix where T=time series length and K=number of covariates (including the intercept)
%	(Y): Vector of length T containing the observed values of the response	
%   -- DLM inputs --	
%	(THETA0): Initial parameters (default: theta=zero) 
%	(S0): Initial parameter estimates (default: 100*eye(dim))
%	(L0): Initial value of forgetting factor, lambda (default: 0.99)
%	(LMAX): Maximum value of forgetting factor (default: 0.999)
%	(LMIN): Minimum value of forgetting factor (default: 0.9)
%	(B): Forgetting factor for observational variance (default: 1)
%   -- ADAM stochastic gradient descent --
%	(A): Initial step-size - IF (A) <=0 Forgetting factor is NOT updated (default: 0.005)
%	(B1): Exponential decay rate for the mean of the gradient (default=0.8)
%	(B2): Exponential decay rate for the variance of the gradient (default=0.8)
%
% Outputs:
%	(Y_PRED): Predicted values for time-series
%	(THETA): Estimated coefficients at each time-step
%	(LAMBDA): Forgetting factor at each time-step
%	(dEdL): Derivative of forecast error w.r.t. Lambda
%	(E_T): 1/2 times squared forecast error
%	(STEP): Step-size computed through ADAM

% Set default values for missing inputs
if nargin <5 | isempty(lambda0), lambda0 = 0.99; end
if nargin <6 | isempty(lambda_max0), lambda_max0 = 0.999; end
if nargin <7 | isempty(lambda_min0), lambda_min = 0.9; end
if nargin <8 | isempty(beta0), beta0=1; end
if nargin <9 | isempty(a), a=0.005; end
if nargin <10 | isempty(b1), b1=0.8; end
if nargin <11 | isempty(b2), b2=0.8; end

[T, dim] = size(X);
if nargin<3 | isempty(theta0), 
	theta0 = zeros(1,dim);
else
	assert(length(theta0)==dim,'Incorrect dimensions for initial theta\n');
end
if nargin<4 | isempty(S0), 
	S0 = 100*eye(dim); 
else
	assert(size(S0,1)==dim & size(S0,2)==dim & issymmetric(S0));
end

assert(length(y)==T,'Length of time series, length(y), incompatible with X');

%%%%%% INITIALISATION
% Forgetting factor
lambda = lambda0*ones(T,1);
lambda_max = lambda_max0;
lambda_min = lambda_min0;
lambdaInv = 1./lambda(1);
beta = beta0;

% Parameter estimates
theta = zeros(T,dim);
theta(1,:) = theta0;

% Predictions
y_pred = zeros(T,1);
%yvar = zeros(T,1);
e_t = zeros(T,1);
%k_t = zeros(dim,T);
V_t = 0;

n = 2;

% covariance matrix of estimated theta (using notation in eDMA JSS paper)
C = S0;

% Initialisation as in eDMA R package: Estimate_par.cpp
R_mat = C;
y_pred(1) = X(1,:) * theta(1,:)';
% one-step-ahead predictive variance
yvar = X(1,:) * R_mat * X(1,:)';

% prediction error
e_t(1) = y(1) - y_pred(1);

%Kalman Gain (vA)
k_t = R_mat*X(1,:)'./yvar;

% Update of theta
theta(1,:) = theta(1,:) + e_t(1) * k_t';  % Eq.(4) in DMA-AF.pdf

% Estimate of process noise variance
V_t = (y(1)^2 + e_t(1)^2/yvar)/n;

% Likelihood
%f = zeros(T,1);

% Recursively estimated derivatives
% dedl - derivative of 1/2(squared forecast error) w.r.t. lambda
dedl = zeros(T,1);
% psi - derivative of theta w.r.t lambda
psi = zeros(dim,1);
% dKdl - derivative of Kalman gain w.r.t. lambda
dKdl = zeros(dim,1);
% Omega - derivative of yvar (one-step-predictive variance) w.r.t. lambda
Omega = 0;
% S - derivative of C (Sigma) w.r.t lambda
S = zeros(dim,dim);
% B - derivative of V (process noise variance) w.r.t lambda
B = 0;

% ADAM related parameters
m = 0;
v = 0;
if nargout >=6, step = zeros(T,1); end
% =============================| end of preliminaries |=========================

% =============================Start the Kalman filter loop   
for t = 2:T
	n = beta*n + 1;
	
	% Update of Variance of estimated theta (Eq. 6 in Korobilis, Eq. 17 in eDMA)
	R_mat = lambdaInv*C;

	% -----------------------Prediction stage
	y_pred(t) = X(t,:) * theta(t-1,:)' ; % predict t given info at t-h
	if isnan(y_pred(t) )
		keyboard
	end

	% Update covariance of predicted y (Eq. 18 in eDMA paper, Eq. 10 in Korobilis)
	yvar = V_t + X(t,:) * R_mat * X(t,:)'; % The forecast variance of each model
	if isnan(yvar)
		fprintf('yvar is NaN');
		keyboard
	end

	% HERE: Update Omega - derivative of yvar w.r.t. lambda
	Omega = B + lambdaInv * X(t,:)*(S - lambdaInv*C)*X(t,:)';
    
	% Prediction error
	e_t(t) = y(t) - y_pred(t); % one-step ahead prediction error, Eq.(3) in DMA-AF.pdf

	% HERE: take gradient of d{e^2}/dl: I need previous step psi
	dedl(t) = -e_t(t) * X(t,:) * psi;

	% Kalman Gain, Eq.(2) in DMA-AF.pdf
	k_t = R_mat*X(t,:)'./yvar;

	% HERE: Update dKdl - derivative of Kalman Gain w.r.t. lambda
	dKdl = lambdaInv*S*X(t,:)'/yvar - (lambdaInv + Omega/yvar) * k_t;

	% HERE: Update B - derivative of V w.r.t lambda
	B = B + (B/n)*(e_t(t)^2/yvar - 1) + ...
		(V_t/n)*(-2*e_t(t)*X(t,:)*psi/yvar - e_t(t)^2 * Omega/yvar^2);

	% New update of process noise variance from eDMA paper
	V_t = V_t + (V_t/n)*(e_t(t)^2/yvar - 1);
   	
	% Update model coefficients
	theta(t,:) = theta(t-1,:) + e_t(t) * k_t';  % Eq.(4) in DMA-AF.pdf 
   
	% HERE: Update psi - derivative of theta w.r.t lambda
	psi = psi + e_t(t)*dKdl - (X(t,:)*psi) * k_t;
    
	% HERE: Update S - derivative of C (Sigma) w.r.t lambda
	S = lambdaInv*S - lambdaInv^2 * C - (k_t* dKdl' + ...
		dKdl*k_t')*yvar - k_t*k_t'*Omega;

	% Update state covariance matrix P_t: Eq. (6) in DMA-AF.pdf
	C = R_mat - k_t * k_t'*yvar; 

	% Update forgetting factor: lambda through ADAM algorithm
	m = b1*m + (1-b1)*dedl(t);
	v = b2*v + (1-b2)*dedl(t)^2;
	if t > 1 & a>0,
		if nargout >=6, 
			step(t) = a*m/((1-b1^t) * (sqrt(v/(1-b2^t)) + 1.e-8));
		end
		lambda(t) = lambda(t-1) - a*m/((1-b1^t) * (sqrt(v/(1-b2^t)) + 1.e-8));
		lambda(t) = max(min(lambda(t), lambda_max),lambda_min);
		lambdaInv = 1./lambda(t);
	end
end
end
