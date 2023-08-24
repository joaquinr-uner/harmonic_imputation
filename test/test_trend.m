
addpath(genpath('..'))
N = 4000;
fs = 4000;
rng(0)
t = 0:1/fs:(N-1)/fs;
%phi = 50*t + 5/(2*pi)*cos(2*pi*t);
phi = -10*t.*log(t+0.1) + 75*t;
%phi = 27*t;
A = ones(1,N);

K = 5;

trend = 0.3*cos(2*4*pi*t);
%trend = zeros(1,N);

x = cos(2*pi*phi) + trend;
subs = {'Const'};
%subs = {'Const','Poly','Bump','Tanh'};
for k=2:K
    haf = sample_haf(t,subs);

    x = x + haf.*cos(2*pi*k*phi);

end

true = A.*x;

true = true - mean(true);

L = round(0.2*N);


true = true';
SNR = Inf;
r = 10^(-SNR/20)*std(x)*randn(N,1);
Lmin = round(L/4);
Lmax = round(L/2);

st1 = round(N/4-0.05*N) + randi(0.1*N);
st2 = round(N/2-0.05*N) + randi(0.1*N);
st3 = round(3*N/4-0.05*N) + randi(0.1*N);

N = length(true);
x = true;

Ll = floor(randfixedsum(3,1,L,Lmin,Lmax)');
ed1 = st1 + Ll(1) - 1;
ed2 = st2 + Ll(2) - 1;
ed3 = st3 + Ll(3) - 1;

x(st1:ed1) = 0;
x(st2:ed2) = 0;
x(st3:ed3) = 0;

[~,fh] = compute_sigma(x);
plot(x)

[sth,Lh] = missing_ints(x,struct('c','x','d',0.01*fs,'t',0));

%params = struct('d',0.5);
params = struct('path','/home/sentey/Dropbox/Github/harmonic_imputation');
imp = impute_arimaf(x,sth,Lh,params);

[AD,phiD,trend_imp] = harm_decomp(imp,struct('with_trend',1));

[s_int,trend_int] = harm_int(AD,phiD,sth,Lh,'pchip',imp,trend_imp);

figure(1)
plot(t,x,'--')
hold on
plot(t,imp,'r')
plot(t,s_int,'k')
hold off


errors_ref = compute_errors(true,x,sth,Lh,{'mae','mse','rmse','sim'});
errors_imp = compute_errors(true,imp,sth,Lh,{'mae','mse','rmse','sim'});
errors_int = compute_errors(true,s_int,sth,Lh,{'mae','mse','rmse','sim'});

