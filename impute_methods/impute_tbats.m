function [s_imp] = impute_tbats(s,st,L,params)
%%
% Missing data imputation by ARIMA model forward forecasting.
% Inputs: 
%         s: signal with missing data.
%         st: initial indexes of missing data intervals.
%         L: lengths of missing data intervals.
%         params: optional parameters struct. Includes:
%                 cycl: number of cycles used to fit ARIMA model
%                 fmax: maximum STFT frequency, as proportion of sampling
%                 frequency.
%                 sigma: window parameter for STFT
% Outputs:
%         s_imp: signal with imputed values.
if nargin<4
    pn = cd;
    T = 0;
    j = 0;
else
    if isfield(params,'pn')
        pn = params.pn;
    else
        pn = cd;
    end
    if isfield(params,'T')
        T = params.T;
    else
        T = 0;
    end
    if isfield(params,'j')
        j = params.j;
    else
        j = 0;
    end
end
s = s(:);

Ni = length(st);

s_imp = s;
if ~isfolder([pn '/temp'])
    mkdir([pn '/temp'])
end
pn_temp = [pn '/temp/'];
intp = 0;
  
for qi=1:Ni

    inti = st(qi):st(qi)+L(qi)-1;

    %sp = s_imp(1:inti(1)-1);
    sp = s_imp(intp(end)+1:inti(1)-1);

    Np = length(sp);
    [sigma,~,b] = compute_sigma(sp);

    [Ff,sFf] = STFT_Gauss(sp,Np,sigma,0.5);

    cf = ridge_ext(Ff,0.1,0.1,10,10);

    [~,phif] = extract_harmonics(Ff,sFf,cf,b,b,1);

    wrapphi = wrapTo2Pi(2*pi*phif);dwphi = diff(wrapphi);
    plocs = find(dwphi<-5);
    pw = diff(plocs);
    if T == 0
        T = floor(median(pw(end-2:end)));
    end
    si = sp(end-sum(pw(end-3:end)):end);
    Ni = length(si);
    sti = st(qi);
    Li = L(qi);
    
    %saveR('temp/MissingSamples.R','si','T','sti','Li','Ni')
    name = tempname(pn_temp);
    saveR([name '.R'],'si','T','sti','Li','Ni','j')
    %CurrentDirectory=strrep(pwd,'\','/');

    evalc(['!/usr/bin/Rscript ' pn '/RunTBATS.R ' name '.R']);

    load([name '.R.mat'],'imp')
        

    s_imp(inti) = imp - mean(imp);
    intp = inti;
    T = 0;
end

%warning('off')
%rmdir(pn_temp,'s')
%warning('on')