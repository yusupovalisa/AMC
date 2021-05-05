function [y_pred, w_star, loss, y_best] = ConfHedge(y, Ypred)
%ConfHedge Fixed Share algorithm 
% [Y_PRED, W, LOSS, Y_BEST] = CONFHEDGE(Y, EXPERTS)
%
% Inputs:
%	(Y): Observed time-series
%	(EXPERTS): Experts' predictions
%
% Outputs:
%	(Y_PRED): Predicted Y at each time-step
%	(W): Weight to each expert
%	(LOSS): Loss of aggregating algorithm at each time-step
%	(Y_BEST): Prediction of expert with highest weight

[T, N] = size(Ypred);
assert(T == length(y));

% Vyugin's parameters
ita = Inf;

% initialise un-normalised weights
%w_t = zeros(T,N);
w_t = 1./N * ones(1,N);

% Normalised weights used for prediction
w_star = zeros(T,N);

% Prediction of aggregating algorithm at time t
%y_t_N = zeros(N,1);

% Cumulative mixloss
Delta_t = 0;

% Loss of every model at each time-step
loss = zeros(T,N);

% Aggregating algorithm prediction
y_pred = zeros(T,1);
% Best model prediction
y_best = zeros(T,1);

for t=1:T
	w_star(t,:) = w_t./sum(w_t);
	
	y_pred(t) = w_star(t,:) * Ypred(t,:)';
	
	% model with highest weight (ensuring that ties are correctly resolved)
	b = max(w_star(t,:));
	id = find(w_star(t,:) == b);
	y_best(t) =  mean(Ypred(t,id));
	
	% Loss at time t for each model
	loss(t,:) = 0.5*(y(t) - Ypred(t,:)).^2;
	
	% Aggregating algorithm loss
	h_t = w_star(t,:) * loss(t,:)';
	
	%Update the weights:
	% Loss Update
	if ~isinf(ita)
		w_mu = w_t .* exp(-ita * loss(t,:));
		% used in mixloss computation
		sumW = sum(w_mu);
		if sumW < eps
			%fprintf('vyugin.m: Sum of Weights effectively zero at iter %i\n', t);
			%keyboard;
		end
		w_mu = w_mu/sumW;
	else
		w_mu = zeros(1,N);
		% Look at p.9 of Follow the Leader If You Can, Hedge If You Must
		mn = min(loss(t,:));
		mnz = find(loss(t,:) == mn);
		w_mu(mnz) = 1;
		sumW = length(mnz);
		w_mu = w_mu/sumW;
	end
	% Preparing weights for next iteration
	% Mixing Scheme: Fixed Share
	alpha = 1/(t+1);
	w_t = alpha/N  + (1-alpha)*w_mu;
	
	if ~isinf(ita)
		m_t = -log(sumW)/ita;
	else
		m_t = 0;
	end
	% see p.9 of Follow the Leader if you can, hedge if you must
	delta_t = max(0,h_t - m_t);
	Delta_t = Delta_t + delta_t;
	ita = max(1, log(N))/Delta_t;
end

% END OF FUNCTION
end

