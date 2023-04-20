
addpath(genpath('..'))
N = 1000;
fs = 1000;
rng(0)
t = 0:1/fs:(N-1)/fs;
%phi = 50*t + 5/(2*pi)*cos(2*pi*t);
phi = 25*t;
A = ones(1,N);

x = cos(2*pi*phi) + 0.5*cos(2*pi*2*phi);

x(N/2-99:N/2+100) = 0;

params = struct('path','/home/sentey/Dropbox/Github/harmonic_imputation');
x_imp = impute_tbats(x,N/2-99,200,params);

figure(1)
plot(t,x,'--')
hold on
plot(t,x_imp)