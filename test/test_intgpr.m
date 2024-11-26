addpath(genpath('C:\Users\Intel\Dropbox\Github'))
addpath(genpath('/home/sentey/Dropbox/Github'))
N = 2000;
fs = 2000;
rng(0)
t = 0:1/fs:(N-1)/fs;
%phi = 50*t + 2.5/(2*pi)*cos(2*pi*t);
phi = 50*t + 5*t.^2;
%phi = 50*t;
A = ones(1,N);

K = 1;

x = cos(2*pi*phi);
subs = {'Const','Bump','Tanh','Trigo'};
%subs = {'Const'};
for k=2:K
    haf = sample_haf(t,subs);

    x = x + haf.*cos(2*pi*k*phi);

end

plot(x)


true = A.*x;

true = true - mean(true);

L = 0.1*N;

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

[sth,Lh] = missing_ints(x,struct('c','x','d',0.01*fs,'t',0));

params_imp = struct();

x_imp = impute(x,sth,Lh,{'tlm'});

[AD,phiD] = harm_decomp(x_imp');
%[AD,phiD] = harm_decomp(x');
tic
kern = 'exponential';
%kern = 'squaredexponential';
%kern = 'matern52';
x_int = harm_int_gpr(AD,phiD,sth,Lh,'gpr',x_imp',zeros(size(x)),kern);
t_intgpr = toc;
figure(1)
plot(x,'k')
hold on
plot(true,'b','LineWidth',1.5)
plot(x_imp,'r')
plot(x_int,'g')
hold off
legend('missing','org','imp','int')

err_ref = compute_errors(true,x,sth,Lh,'mae');
err_imp = compute_errors(true,x_imp,sth,Lh,'mae');
err_igpr = compute_errors(true,x_int,sth,Lh,'mae');