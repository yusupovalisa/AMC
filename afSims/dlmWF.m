function [F, y, theta_hat, theta, epsilon] = dlmWF(lambda,T0,p0,V,C0,Seed)
% Simulate data using DLM with forgetting/ discounting 
% Linear system representation
% theta_t = theta_t-1 + ita_t
% y_t = F_t'*theta_t + e_t
% with:
% ita_t = N(0,W), and e_t = N(0,V),
% where
% F - observation/ covariates matrix
% W - state noise covariance
% V - observation noise covariance

% State transition noise
if nargin>5 & ~isempty(Seed),
	rng(seed);
end

beta = 1;
% Time periods
T = T0;
% Number of time-varying regression coefficients
p = p0; 

% X matrix (using notation in eDMA JSS paper)
F = randn(T, p);
y = zeros(T,1);

% Initial value of theta covariance matrix
C = C0;

% Measurement noise at every time-step
epsilon = sqrt(V)*randn(T,1);

% Fixed value of lambda
lambdaInv = 1/lambda;

% The state noise variance
W = (1-lambda)*lambdaInv*C;
% %The state noise.
ita = chol(W) * randn(p0,1);

% Initial condition of the state
% Fix the values of theta
theta = zeros(p,T);
theta_hat = zeros(p,T);
theta_0 = zeros(p,1);%chol(0.1*C)*randn(p,1);%zeros(p,1); %chol(0.1*C)*randn(p,1); %zeros(p,1);

% Initial state and initial value of y
theta(:,1) = theta_0 + ita;
y(1) = F(1,:)*theta(:,1) + epsilon(1);

% Update equations
% Predictive variance of y
n=2;
Q = F(1,:)*C*F(1,:)';
% Prediction error
e(1) = y(1) - F(1,:)*theta_hat(:,1);
% Kalman gain
k(:,1)= C*F(1,:)'./Q;

% Covariance matrix of theta
%C= lambdaInv*C - k(:,1)*F(1,:)*lambdaInv*C;
% Theta update
theta_hat(:,1) = theta_hat(:,1) + e(1)*k(:,1);
V_t = (y(1)^2 + e(1)^2/Q)/n;

for t = 2:T,
	n = beta*n + 1;
	% State noise covariance
	W = (1-lambda)*lambdaInv*C;
	% State noise.
	ita = chol(W) * randn(p0,1);
	
	% State transition equation
	theta(:,t) = theta(:,t-1) + ita;

	% Measurement equation
	y(t) = F(t,:)*theta(:,t) + epsilon(t);

	%Prediction error
	e(t) = y(t) - F(t,:)*theta_hat(:,t-1);

	%Updating equations
	% Predictive variance of y
	Q = V_t + lambdaInv* F(t,:)* C * F(t,:)';

	V_t = V_t + (V_t/n)*(e(t)^2/Q - 1);

	% Kalman gain
	k(:,t) = (lambdaInv*C*F(t,:)')./Q;
	% Covariance matrix of theta
	%C = lambdaInv*C - k(:,t)*F(t,:)*lambdaInv*C;
	C = lambdaInv*C - k(:,t)*k(:,t)'*Q;

	% Theta update
	theta_hat(:,t) = theta_hat(:,t-1)+k(:,t)*e(t);   
end    

theta = theta';
theta_hat = theta_hat';

% plot(y)
% hold on
% plot(theta(1,:))
% plot(theta_hat(1,:))
% hold off
