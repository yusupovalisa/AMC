function [y_pred, w, lgL] = opp(y, Ypred, Fpred)

[T, N] = size(Fpred);
assert(T == length(y));

% Aggregating algorithm prediction
y_pred = zeros(T,1);
% Best model prediction
lgL = zeros(T,1);

w = zeros(T,N);
w(1,:) = ones(1,N)./N;

y_pred(1) = w(1,:) * Ypred(1,:)';
lgL(1) = -logL(y(1), w(1,:), Fpred);

for t=2:T,
	%keyboard
	options = optimoptions('fmincon','SpecifyObjectiveGradient',true, ...
		'Display','off');%,'Algorithm','sqp');
	fun = @(w)(logL(y(1:t-1),w, Fpred));
	[x,fval,flag,output] = fmincon(fun,w(t-1,:),[],[],ones(1,N),1, ...
		zeros(1,N),ones(1,N),[],options);

	if flag < 0,
		fprintf('HMweights: failed to identify minimiser\n');
		%keyboard;
	end
	w(t,:) = x;
	lgL(t) = -t*fval + sum(lgL(1:t-1));
	y_pred(t) = w(t,:) * Ypred(t,:)';
end
end

function [F,G] = logL(y1, w1, F1pred)
% returns Negative log-Likelihood for optimisation purposes
N = length(w1);
T = length(y1);

F = 0;
G = zeros(N,1);
g = zeros(N,1);
for t=1:T,
	%f = 0;
	for i=1:N,
		g(i) = F1pred{t,i}(y1(t));
		%f = f + w1(t,i)*g(i);
	end
	f = w1*g;
	%g = g./f;
	G = G + g./f;
	F = F + log(f);
end
F = -F/T;
G = -G/T;
end
