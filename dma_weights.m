function [y_pred, w, y_best] = dma_weights(y, Mpred, Spred, a, C)

if nargin < 5, C = 0; end

T = length(y);
K = size(Mpred,2);
L = zeros(T, K);
w = zeros(T, K);
w(1,:) = ones(1,K)/K;
y_pred = zeros(T,1);
y_best = zeros(T,1);
for t=1:T-1,
	y_pred(t) = w(t,:) * Mpred(t,:)';

	id = find(w(t,:) == max(w(t,:)));
	y_best(t) =  mean(Mpred(t,id));

	for k=1:K,
		L(t,k) = normpdf(y(t), Mpred(t,k), Spred(t,k));
	end
	S = sum(w(t,:).^a + C);
	w(t+1,:) = L(t,:) .* ((w(t,:).^a + C)/S);
	w(t+1,:) = w(t+1,:)./sum(w(t+1,:));
end
