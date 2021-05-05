addpath('afSims/')

%%%%%%%%%%%% ADAPTIVE FORGETTING/ DISCOUNTING SIMULATIONS
Seed = 1989;
% ABRUPT CHANGE
rng(Seed);
T=1000;
iRep = 100;
theta0 = [3;2;1;-1;-2];
d = 5;

y_pred = zeros(T,iRep);
e_t = y_pred;
lambda = y_pred;
thetaHat = zeros(T,iRep*5);
ctime = 0;
for r = 1:iRep,
	%fprintf('Iteration %i\n', r);
	[y, X, theta] = dlm_abrupt(T,5,1.0,theta0);

	tic;
	[y_pred(:,r), thetaHat(:,[r:100:500]), lambda(:,r), e_t(:,r)] = afdlm(X, y, zeros(d,1), 100*eye(d), ...
		0.99, 0.999, 0.9, 1, 0.003, 0.8, 0.8);
	ctime = ctime + toc;
end
fprintf('Average time %1.3f\n', mean(ctime));
% 1/2 ( squared forecast error )
e_t = 0.5*e_t.^2;

results = [mean(y_pred,2), quantile(y_pred,[0.5,0.25,0.75,0.025,0.975],2),...
	mean(e_t,2), quantile(e_t,[0.5,0.25,0.75,0.025,0.975],2),...
	mean(thetaHat(:,1:iRep),2), quantile(thetaHat(:,1:iRep),[0.5,0.25,0.75,0.025,0.975],2),...
	mean(thetaHat(:,101:2*iRep),2), quantile(thetaHat(:,101:2*iRep),[0.5,0.25,0.75,0.025,0.975],2),...
	mean(thetaHat(:,201:3*iRep),2), quantile(thetaHat(:,201:3*iRep),[0.5,0.25,0.75,0.025,0.975],2),...
	mean(thetaHat(:,301:4*iRep),2), quantile(thetaHat(:,301:4*iRep),[0.5,0.25,0.75,0.025,0.975],2),...
	mean(thetaHat(:,401:5*iRep),2), quantile(thetaHat(:,401:5*iRep),[0.5,0.25,0.75,0.025,0.975],2),... 
	mean(lambda,2), quantile(lambda,[0.5,0.25,0.75,0.025,0.975],2)];

Tb = array2table([[1:T]', theta', results], 'VariableNames', {'time','th1','th2','th3','th4','th5',...
	'muYhat', 'medYhat','q1Yhat','q3Yhat','q025Yhat','q975Yhat',...
	'muEr', 'medEr','q1Er','q3Er','q025Er','q975Er',...
	'muth1', 'medth1','q1th1','q3th1','q025th1','q975th1',...
	'muth2', 'medth2','q1th2','q3th2','q025th2','q975th2',...
	'muth3', 'medth3','q1th3','q3th3','q025th3','q975th3',...
	'muth4', 'medth4','q1th4','q3th4','q025th4','q975th4',...
	'muth5', 'medth5','q1th5','q3th5','q025th5','q975th5', ...
	'muL', 'medL','q1L','q3L','q025L','q975L'});
writetable(Tb,'afSims/AFDLM_sim_abrupt.csv');

% GRADUAL DRIFT
rng(Seed);
% initialise file
Tb = table([],[],[],[],[],[],[],[], 'VariableNames', {'time', 'lambda', 'muL', 'medL','q1L','q3L','q025L','q975L'});
writetable(Tb, 'afSims/AFDLM_sim_gradual.csv');
for TrueL = [0.99, 0.97, 0.95], 
	for r = 1:iRep,
		[X, y] = dlmWF(TrueL, T, d, 1.0, 0.1*eye(d));

		[y_pred, thetaHat, lambda(:,r)] = afdlm(X, y, zeros(d,1), 0.1*eye(d), ...
			0.99, 0.999, 0.9, 1, 0.003, 0.8, 0.8);
	end
	results = [[1:T]', TrueL*ones(T,1), mean(lambda,2), quantile(lambda,[0.5,0.25,0.75,0.025,0.975],2)];
	Tb = array2table(results, 'VariableNames', {'time','lambda','muL','medL','q1L','q3L','q025L','q975L'});
	writetable(Tb, 'afSims/AFDLM_sim_gradual.csv', 'WriteMode','Append');
end


 STATIC 
rng(Seed);
T=1000;
iRep = 100;
theta0 = [3;2;1;-1;-2];
d = 5;
y_pred = zeros(T,iRep);
lambda = y_pred;
thetaHat = zeros(T,iRep*5);
for r = 1:iRep,
	%fprintf('Iteration %i\n', r);
	X = randn(T, d);
	y = X*theta0 + randn(T,1);
	[y_pred, thetaHat(:,[r:100:500]), lambda(:,r)] = afdlm(X, y, zeros(d,1), 100*eye(d), ...
		0.99, 0.999, 0.9, 1, 0.003, 0.8, 0.8);
end
results = [ [1:T]', repmat(theta0',T,1), ...
	mean(thetaHat(:,1:iRep),2), quantile(thetaHat(:,1:iRep),[0.5,0.25,0.75,0.025,0.975],2),...
	mean(thetaHat(:,101:2*iRep),2), quantile(thetaHat(:,101:2*iRep),[0.5,0.25,0.75,0.025,0.975],2),...
	mean(thetaHat(:,201:3*iRep),2), quantile(thetaHat(:,201:3*iRep),[0.5,0.25,0.75,0.025,0.975],2),...
	mean(thetaHat(:,301:4*iRep),2), quantile(thetaHat(:,301:4*iRep),[0.5,0.25,0.75,0.025,0.975],2),...
	mean(thetaHat(:,401:5*iRep),2), quantile(thetaHat(:,401:5*iRep),[0.5,0.25,0.75,0.025,0.975],2),... 
	mean(lambda,2), quantile(lambda,[0.5,0.25,0.75,0.025,0.975],2)];

Tb = array2table(results, 'VariableNames', {'time','th1','th2','th3','th4','th5',...%	'muth1', 'medth1','q1th1','q3th1','q025th1','q975th1',...
	'muth2', 'medth2','q1th2','q3th2','q025th2','q975th2',...
	'muth3', 'medth3','q1th3','q3th3','q025th3','q975th3',...
	'muth4', 'medth4','q1th4','q3th4','q025th4','q975th4',...
	'muth5', 'medth5','q1th5','q3th5','q025th5','q975th5', ...
	'muL', 'medL','q1L','q3L','q025L','q975L'});
writetable(Tb,'afSims/AFDLM_sim_static.csv');


%%%%%%%%%%%% MODEL COMBINATION SIMULATIONS
addpath('mcSims/')
%% Persistence
%iRep = 100;
%T = 300;
%Weights = zeros(T, iRep);
%
%Vnames = {'phi','time', 'medW1', 'q1W1', 'q3W1', 'medW2','q1W2','q3W3', 'medW3','q1W3','q3W2'};
%Tb = table([],[],[],[],[],[],[],[],[],[],[], 'VariableNames', Vnames);
%writetable(Tb, 'mcSims/ConfHedge_persistence.csv');
%rng(Seed);
%for phi=[0.01:0.01:0.99],
%	for r = 1:iRep,
%		[y, Ypred] = persistence(T, phi);
%		[y_v, Weights(:,r:iRep:3*iRep)] = ConfHedge(y, Ypred);
%	end
%	results = [phi*ones(T,1), [1:T]', quantile(Weights(:,1:iRep),[0.5,0.25,0.75],2), ...
%		quantile(Weights(:,(iRep+1):2*iRep),[0.5,0.25,0.75],2), ...
%		quantile(Weights(:,(2*iRep+1):3*iRep),[0.5,0.25,0.75],2)];
%	Tb = array2table(results, 'VariableNames', Vnames);
%	writetable(Tb, 'mcSims/ConfHedge_persistence.csv', 'WriteMode','Append');
%end

% Model combination under abrupt change:
Vnames = {'alg', 'time', 'error', 'w1', 'q1w1','q3w1', 'w2','q1w2','q3w2','w3','q1w3','q3w3'};
Tb = table([],[],[],[],[],[],[],[],[],[],[],[], 'VariableNames',Vnames);
writetable(Tb, 'mcSims/abrupt_change.csv');
T=300;
% BMA/ DMA
er = zeros(T,iRep*3);
wDMA = zeros(T,iRep*3*3);
rng(Seed);
for r = 1:iRep,
	[y, Ypred, Fpred, Spred] = mcAbrupt();
	c = 0;
	for lambda = [1, 0.99, 0.95],
		cols = (c*3*iRep+r):iRep:(c+1)*3*iRep;	
		[y_DMA, wDMA(:, cols)] = dma_weights(y, Ypred, Spred, lambda, 0);
		er(:,(c*iRep+r)) = 0.5*(y-y_DMA).^2;
		c = c+1;
	end
end
Names = {'BMA','DMA99','DMA95'};
for c = 0:2,
	cols = (c*3*iRep+1):(c*3+1)*iRep;
	results=[[1:T]', median(er(:,(c*iRep+1):(c+1)*iRep),2), ...
		quantile(wDMA(:,cols),[0.5,0.25,0.75],2), ...
		quantile(wDMA(:,cols+iRep),[0.5,0.25,0.75],2), ...
		quantile(wDMA(:,cols+2*iRep),[0.5,0.25,0.75],2)];

	Tab = array2table(results, 'VariableNames', {Vnames{2:end}});
	alg = cellstr(repmat(Names{c+1},T,1));
	Tab = addvars(Tab, alg ,'Before','time');
	writetable(Tab, 'mcSims/abrupt_change.csv','WriteMode','Append');
end

% CONFHEDGE
er = zeros(T,iRep);
w = zeros(T,3*iRep);
rng(Seed);
for r = 1:iRep,
	fprintf('OPP iteration %i\n', r)
	[y, Ypred, Fpred, Spred] = mcAbrupt();
	[y_v, w(:,[r:iRep:3*iRep])] = ConfHedge(y, Ypred);
	er(:,r) = 0.5*(y - y_v).^2;
end
Tab = array2table([[1:T]', median(er,2), ...
	quantile(w(:,1:iRep),[0.5,0.25,0.75],2), ...
	quantile(w(:,(1+iRep):2*iRep),[0.5,0.25,0.75],2), ...
	quantile(w(:,(1+2*iRep):3*iRep),[0.5,0.25,0.75],2)], 'VariableNames',{Vnames{2:end}});
alg = cellstr(repmat('ConfHedge',T,1));
Tab = addvars(Tab, alg ,'Before','time');
writetable(Tab, 'mcSims/abrupt_change.csv','WriteMode','Append');

% Optimal Prediction Pool
er = zeros(T,iRep);
w = zeros(T,3*iRep);
rng(Seed);
for r = 1:iRep,
	[y, Ypred, Fpred, Spred] = mcAbrupt();
	[y_opp, w(:,[r:iRep:3*iRep])] = opp(y, Ypred, Fpred);
	er(:,r) = 0.5*(y - y_opp).^2;
end
Tab = array2table([[1:T]', median(er,2), ...
	quantile(w(:,1:iRep), [0.5,0.25,0.75],2), ...
	quantile(w(:,(1+iRep):2*iRep), [0.5,0.25,0.75],2), ...
	quantile(w(:,(1+2*iRep):3*iRep), [0.5,0.25,0.75],2)], 'VariableNames',{Vnames{2:end}});
alg = cellstr(repmat('OPP',T,1));
Tab = addvars(Tab, alg ,'Before','time');
writetable(Tab, 'mcSims/abrupt_change.csv','WriteMode','Append');

% Reproduction of DeCo and BPS results requires the corresponding packages, which can be found at:
% DECO: https://www.jstatsoft.org/article/view/v068i03/
% BPS: https://www2.stat.duke.edu/~mw/mwsoftware/BPS/index.html
