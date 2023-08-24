addpath(genpath('C:\Users\Intel\Dropbox\Github'))
fs = 4000;
N = 4000;
rng(0)
t = 0:1/fs:(N-1)/fs;
f0 = 50;
%phi = f0*t + 5/(2*pi)*cos(2*pi*t);
phi = f0*t;
A = ones(1,N);

K = 2;

x = cos(2*pi*phi);
subs = {'Const'};
%subs = {'Const','Poly','Bump','Tanh'};
%e = [1 1.995 3.005 4.002 4.999];
e = 1:K;
for k=2:K
    haf = sample_haf(t,subs);

    x = x + haf.*cos(2*pi*e(k)*phi);

end
Mdl = arima('Constant',0.5,'AR',{0.7 0.25},'Variance',.1);

%x = simulate(Mdl,N)';
plot(x)

trend = zeros(1,N);

%true = A.*x+0.1*std(x)*randn(1,N);
true = A.*x;

true = true + trend;

true = true - mean(true);

L = round(0.1*N);

true = true';
SNR = Inf;
r = 10^(-SNR/20)*std(x)*randn(N,1);
Lmin = round(L/4);
Lmax = round(L/2);

st = 0.5*N;
%st1 = round(N/4-0.05*N) + randi(0.1*N);
%st2 = round(N/2-0.05*N) + randi(0.1*N);
%st3 = round(3*N/4-0.05*N) + randi(0.1*N);

N = length(true);
x = true;

ed = st + L - 1;
%Ll = floor(randfixedsum(3,1,L,Lmin,Lmax)');
%ed1 = st1 + Ll(1) - 1;
%ed2 = st2 + Ll(2) - 1;
%ed3 = st3 + Ll(3) - 1;

x(st:ed) = 0;
%x(st1:ed1) = 0;
%x(st2:ed2) = 0;
%x(st3:ed3) = 0;

[~,fh] = compute_sigma(x);
plot(x)

[sth,Lh] = missing_ints(x,struct('c','x','d',0.01*fs,'t',0));

%params = struct('pn','/home/sentey/Dropbox/Github/harmonic_imputation/impute_methods/aux-functs');
params = struct();
%params = struct('d',0.5);
tic
imp = impute_lsw(x,sth,Lh,params);
t_imp = toc;

errors_ref = compute_errors(true,x,sth,Lh,{'mae','mse','rmse','sim'});
errors = compute_errors(true,imp,sth,Lh,{'mae','mse','rmse','sim'});
figure(1)
plot(t,true,'r-')
hold on
plot(t,x,'k--')
plot(t,imp,'b-')
legend('true','missing','imputed')
hold off