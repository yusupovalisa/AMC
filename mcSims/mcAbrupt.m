function [y, Ypred, Fpred, Spred, mX] = mcAbrupt()

% Unconditional means and stds for each time series
m1 = 1; s1 = 2.;
m2 = 0.; s2 = 0.8;
m3 = -2.5; s3 = 0.3;
mu = [m1,m2,m3];
sd = [s1,s2,s3];
% slopes; intercepts and noise variance for each AR model based on formulae:
% Unconditional mean = b0/(1-b1); Unconditional variance = sigma^2/(1-b1^2)
p1 = 0.25  ; c1 = m1*(1-p1) ; v1 = s1^2 *(1-p1^2) ; 
p2 = 0.375 ; c2 = m2*(1-p2) ; v2 = s2^2 *(1-p2^2) ; 
p3 = 0.56  ; c3 = m3*(1-p3) ; v3 = s3^2 *(1-p3^2) ; 
c = [c1,c2,c3];
v = [v1,v2,v3];
p = [p1,p2,p3];

T=300;
Y = zeros(T+1,3);
for t=2:T+1,
	Y(t,:) = c + Y(t-1,:).*p + randn(1,3).*sqrt(v);
end

Ypred = zeros(T,3);
Spred = [sqrt(v1)*ones(T,1), sqrt(v2)*ones(T,1), sqrt(v3)*ones(T,1)];
Fpred = cell(T,3);
if nargout == 5,
	mX = zeros(T,1000,3);
end

Y = zeros(T+1,3);
for t=2:T+1,
        Ypred(t-1,:) = c + Y(t-1,:).*p;
        Y(t,:) = Ypred(t-1,:) + randn(1,3).*sqrt(v);
        for m=1:3,
                Fpred{t-1,m} = @(x)(normpdf(x, Ypred(t-1,m), Spred(m,1)));
		% Relevant for DeCo
		if nargout == 5,
			mX(t-1,:,m) = randn(1,1000)*Spred(m,1) + Ypred(t-1,m);
		end
        end
end

% observed time series
y = [Y(2:101,1); Y(102:201,2); Y(202:301,3)];
