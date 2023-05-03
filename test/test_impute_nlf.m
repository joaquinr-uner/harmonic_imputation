
addpath(genpath('..'))
N = 4000;
fs = 2000;
rng(0)
t = 0:1/fs:(N-1)/fs;
phi = 50*t + 5/(2*pi)*cos(2*pi*t);
%phi = 100*t;
A = ones(1,N);

K = 2;

x = cos(2*pi*phi);
%subs = {'Const','Bump','Tanh','Trigo'};
subs = {'Const'};
for k=2:K
    haf = sample_haf(t,subs);

    x = x + haf.*cos(2*pi*k*phi);

end

plot(x)


true = A.*x;

true = true - mean(true);

L = 800;

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

%x = true + 10^(-10/20)*std(x)*randn(1,N)';

Ll = floor(randfixedsum(3,1,L,Lmin,Lmax)');
ed1 = st1 + Ll(1) - 1;
ed2 = st2 + Ll(2) - 1;
ed3 = st3 + Ll(3) - 1;

x(st1:ed1) = 0;
x(st2:ed2) = 0;
x(st3:ed3) = 0;

plot(x)

[sth,Lh] = missing_ints(x,0.01*fs,0);

[~,fh] = compute_sigma(true);
Th = N/fh;
%M = 100;
%K = floor(2.5*M);
%M = floor(0.9*Th);
%K = ceil(3*Th);
%M = floor(3*Th);
%K = floor(2.5*M);
M = 0;K = 0;
params = struct('M',M,'K',K);
imp_lse = impute_lse(x,sth,Lh,params);
imp_dmd = impute_dmd(x,sth,Lh,params);
imp_gpr = impute_gpr(x,sth,Lh,params);

errors_ref = compute_errors(true,x,sth,Lh,{'mae','mse','rmse','sim'});
errors_lse = compute_errors(true,imp_lse,sth,Lh,{'mae','mse','rmse','sim'});
errors_dmd = compute_errors(true,imp_dmd,sth,Lh,{'mae','mse','rmse','sim'});
errors_gpr = compute_errors(true,imp_gpr,sth,Lh,{'mae','mse','rmse','sim'});

figure(1)
subplot(221)
plot(t,true,'r-')
hold on
plot(t,x,'k--')
plot(t,imp_lse,'b-')
legend('true','missing','imputed LSE')
hold off
subplot(222)
plot(t,true,'r-')
hold on
plot(t,x,'k--')
plot(t,imp_dmd,'b-')
legend('true','missing','imputed DMD')
hold off
subplot(223)
plot(t,true,'r-')
hold on
plot(t,x,'k--')
plot(t,imp_gpr,'b-')
hold off
legend('true','missing','imputed GPR')
hold off