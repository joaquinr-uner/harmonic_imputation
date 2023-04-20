
addpath(genpath('..'))
N = 4000;
fs = 4000;
rng(0)
t = 0:1/fs:(N-1)/fs;
phi = 50*t + 5/(2*pi)*cos(2*pi*t);
%phi = 50*t;
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

L = 0.2*N;
Ni = 3;
[st, ed, Ll] = rand_missing_ints(L,N,Ni);

true = true';
SNR = Inf;
r = 10^(-SNR/20)*std(x)*randn(N,1);

for i=1:Ni
    x(st(i):ed(i)) = 0;
end

[~,fh] = compute_sigma(x);
fh = fh*fs/N;
plot(x)

[sth,Lh] = missing_ints(x,0.01*fs,0);

params = struct('D',K,'fh',fh,'fs',fs);
imp = impute_ddtfa(x,sth,Lh);

errors_ref = compute_errors(true,x,sth,Lh,{'mae','mse','rmse','sim'});
errors = compute_errors(true,imp,sth,Lh,{'mae','mse','rmse','sim'});
figure(1)
plot(t,true,'r-')
hold on
plot(t,x,'k--')
plot(t,imp,'b-')
legend('true','missing','imputed')
hold off