
N = 12000;
fs = 4000;
rng(0)
t = 0:1/fs:(N-1)/fs;
phi = 50*t + 5/(2*pi)*cos(2*pi*t);
%phi = 50*t;
A = ones(1,N);

K = 1;

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

L = 0.05*N;

true = true';
SNR = Inf;
r = 10^(-SNR/20)*std(x)*randn(N,1);
Lmin = round(L/4);
Lmax = round(L/2);

st = 0.5*N;
st1 = round(N/4-0.05*N) + randi(0.1*N);
st2 = round(N/2-0.05*N) + randi(0.1*N);
st3 = round(3*N/4-0.05*N) + randi(0.1*N);

N = length(true);
x = true;

%x = true + 10^(-10/20)*std(x)*randn(1,N)';

Ll = floor(randfixedsum(3,1,L,Lmin,Lmax)');
ed = st + L - 1;
ed1 = st1 + Ll(1) - 1;
ed2 = st2 + Ll(2) - 1;
ed3 = st3 + Ll(3) - 1;

x(st:ed) = 0;
%x(st1:ed1) = 0;
%x(st2:ed2) = 0;
%x(st3:ed3) = 0;

plot(x)

[sth,Lh] = missing_ints(x,struct('c','x','d',0.01*fs,'t',0));

[~,fh] = compute_sigma(true);
Th = N/fh;
%M = 100;
%K = floor(2.5*M);
%M = floor(3*Th);
%K = floor(0.9*Th);
%M = floor(3*Th);
%K = floor(2.5*M);
M = 0;K = 0;
for i=1:10
    sigma = 10^i;
    params = struct('M',M,'K',K,'sigma',sigma);
    imp_edmd = impute_edmd(x,sth,Lh,params);

    errors_ref = compute_errors(true,x,sth,Lh,{'mae','mse','rmse','sim'});
    errors_edmd = compute_errors(true,imp_edmd,sth,Lh,{'mae','mse','rmse','sim'});

    figure(1)
    plot(t,true,'r-')
    hold on
    plot(t,x,'k--')
    plot(t,imp_edmd,'b-')
    legend('true','missing','imputed EDMD')
    title(num2str(sigma))
    hold off
    pause()
end