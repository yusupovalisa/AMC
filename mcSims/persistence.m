function [y, Ypred] = persistence(T, phi)
%This particular experiment was first proposed by Billio et al. (2013) Journal of Econometrics

b0 = [0.1, 0.125, 0.5];
b1 = [phi, 0.5, 0.2];
sd = [0.05, 0.05, 0.05];

Ypred = zeros(T,3);
Y = [[0.25, 0.25, 0.25]; zeros(T,3)];
for t=2:T+1,
        Ypred(t-1,:) = b0 + b1.*Y(t-1,:);
	Y(t,:) =  Ypred(t-1,:) + sd.*randn(1,3);
end
y = Y(2:end,1);
