
addpath(genpath('..'))
addpath(genpath('/home/sentey/Dropbox/Github/MissingData_Real'))
N = 4000;
fs = 4000;
rng(0)
t = 0:1/fs:(N-1)/fs;
phi = 50*t + 5/(2*pi)*cos(2*pi*t);
%phi = 100*t;
A = ones(1,N);

K = 2;

x = cos(2*pi*phi);
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

params = struct();
impf = impute_arimaf(x,sth,Lh,params);
impb = impute_arimab(x,sth,Lh,params);

errors_ref = compute_errors(true,x,sth,Lh,{'mae','mse','rmse','sim'});
errors_arimaf = compute_errors(true,impf,sth,Lh,{'mae','mse','rmse','sim'});
errors_arimab = compute_errors(true,impb,sth,Lh,{'mae','mse','rmse','sim'});

figure(1)
subplot(211)
plot(t,true,'r-')
hold on
plot(t,x,'k--')
plot(t,impf,'b-')
hold off
legend('true','missing','imputed (fw)')
subplot(212)
plot(t,true,'r-')
hold on
plot(t,x,'k--')
plot(t,impb,'g-')
hold off
legend('true','missing','imputed (bw)')
hold off