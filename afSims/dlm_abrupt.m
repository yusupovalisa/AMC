function [y, X, theta] = dlm_abrupt(T,d,V,theta0, Seed)

if nargin>4 & ~isempty(Seed)
	rng(Seed);
end


y = zeros(T,1);

% Values of covariates at all time-points
X = randn(T, d);
% noise
epsilon = sqrt(V)*randn(T,1); 

theta = zeros(d,T);
theta(:,1) = theta0;

y(1) = X(1,:)*theta(:,1) + epsilon(1);
for t = 2:T
    % abrupt change point
    if t == 100,
	    theta(:,t) = theta(:,t-1)*0.5;
    elseif t==400, 
	    theta(:,t) = theta(:,t-1)*1.4;
    elseif t==700,    
	    theta(:,t) = theta(:,t-1)*0.7;
    else
	    theta(:,t) = theta(:,t-1);
    end

    % Measurement equation
    y(t) = X(t,:)*theta(:,t) + epsilon(t);
end    
