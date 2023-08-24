addpath(genpath('/home/sentey/Dropbox/Github/MissingData_Synth'))
addpath(genpath('..'))
N = 8000;
fs = 8000;
rng(0)
t = 0:1/fs:(N-1)/fs;
f0 = 50;
phi = f0*t + (f0/10)/(2*pi)*cos(2*pi*t);
%phi = f0*t;
%A = ones(1,N);
A = 2-cos(2*pi*t);

K = 5;

x = cos(2*pi*phi);
%subs = {'Const'};
subs = {'Const','Poly','Bump','Tanh'};
for k=2:K
    haf = sample_haf(t,subs);

    x = x + haf.*cos(2*pi*k*phi);

end


trend = zeros(1,N);

%true = x;
true = A.*x;

true = true + trend;

true = true - mean(true);

true = true';
SNR = Inf;
r = 10^(-SNR/20)*std(x)*randn(N,1);

M = 0.65*max(true);

x = true;
x(x>M) = M;

figure(1)
plot(t,x)

[sth,Lh] = missing_ints(x,struct('c','y','d',2));

params_decomp = struct();

params = struct('K',2,'M',5);
imp = impute_gpr(x,sth,Lh,params);
[ADms,phiDms] = harm_decomp(x,params_decomp);
[AD,phiD] = harm_decomp(x,params_decomp);
AD = AD(1:K,:);
phiD = phiD(1:K,:);
int = harm_int(AD,phiD,sth,Lh,'spline',x);
intov = int_oversat(x,sth,Lh,'spline');

errors_ref = compute_errors(true,x,sth,Lh,{'mae','mse','rmse','sim'});
errors_imp = compute_errors(true,imp,sth,Lh,{'mae','mse','rmse','sim'});
errors_int = compute_errors(true,int,sth,Lh,{'mae','mse','rmse','sim'});
errors_intov = compute_errors(true,intov,sth,Lh,{'mae','mse','rmse','sim'});
figure(2)
plot(t,true,'r-')
hold on
plot(t,x,'k--')
plot(t,imp,'b-')
plot(t,int,'g--')
plot(t,intov,'c--')
legend('true','missing','imputed','interpolated')
hold off
