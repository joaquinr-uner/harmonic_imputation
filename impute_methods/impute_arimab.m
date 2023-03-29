function [s_imp] = impute_arimab(s,st,L,params)
%%
% Missing data imputation by ARIMA model backward forecasting.
% Inputs: 
%         s: signal with missing data.
%         st: initial indexes of missing data intervals.
%         L: lengths of missing data intervals.
%         params: optional parameters struct. Includes:
%                 cycl: number of cycles used to fit ARIMA model
%                 b: harmonic reconstruction half-width
%                 sigma: STFT gaussian window length
%                 fmax: maximum STFT frequency, as proportion of sampling
%                 frequency.
% Outputs:
%         s_imp: signal with imputed values.
if nargin<4
    cycl = 3;
    sigma = compute_sigma(s,2);
    b = round(sqrt(log(80)*sigma/pi^2)*length(s)) + 1;
    fmax = 0.5;
else
    if isfield(params,'cycl')
        cycl = params.cycl;
    else
        cycl = 3;
    end
    if isfield(params,'sigma')
        sigma = params.sigma;
    else
        sigma = compute_sigma(x,fs,2);
    end
    if isfield(params,'b')
        b = params.b;
    else
        b = round(sqrt(log(80)*sigma/pi^2)*length(s)) + 1;
    end
    if isfield(params,'fmax')
        fmax = params.fmax;
    else
        fmax = 0.5;
    end
    if isfield(params,'redun')
        redun = params.redun;
    else
        redun = 1;
    end
end
s = s(:);

N = length(s);
Ni = length(st);

s_imp = s;

for qi=Ni:-1:1

    inti = st(qi)+L(qi)+1:N;

    %sp = s_imp(1:inti(1)-1);
    sp = s_imp(intp(end)+1:inti(1)-1);
    
    [Ff,sFf] = STFT_Gauss(sp,length(sp)*redun,sigma,fmax);

    cf = ridge_ext(Ff,0.1,0.1,10,10);

    [~,phif] = extract_harmonics(Ff,sFf,cf,b,b,1);

    %phif = phif/redun;
    Np = L(qi);

    ext_arima = extendSig(sp,phif,cycl,Np,'bw');

    s_imp(st(qi)+1:end) = ext_arima;
    intp = inti;
end
