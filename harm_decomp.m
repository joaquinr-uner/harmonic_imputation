function [AD,phiD,TFdata] = harm_decomp(x,params)
%%
% Compute the harmonic decomposition of signal x.
% Inputs: 
%         x: signal under analysis
%         fs: sampling frequency of x
%         params: optional parameters struct. Includes:
%                 'sigma': window parameter sigma.
%                 'b': half-width for harmonic reconstruction
%                 'fmax': maximum frequency for the computation of the
%                 STFT, a proportion of fs
%                 'r_max': maximum admissible number of harmonics
%                 'Criteria': string or cell with criteria for automatic
%                 number of harmonic estimation.
%                 'crit_params': Additional parameters for trigonometric
%                 regression criteria:
%                       When 'Criteria' == 'Wang': vector of penalization
%                                                  parameters c.
%                       When 'Criteria' == 'Kavalieris': maximum order H.
%                       When 'Criteria' == 'Rl': noise variance estimate
%                       sigma.

if nargin<3
    sigma = compute_sigma(x,1);
    %b = round(sqrt(log(80)*sigma/pi^2)*length(x)) + 1;
    b = round(3/pi*sqrt(sigma/2)*length(x));
    fmax = 0.5;
    r_max = 50;
    crit = {'Wang'};
    crit_params = struct('c',[4,6,8,10,12]);
else
    if isfield(params,'sigma')
        sigma = params.sigma;
    else
        sigma = compute_sigma(x,1);
    end
    if isfield(params,'b')
        b = params.b;
    else
        %b = round(sqrt(log(80)*sigma/pi^2)*length(x)) + 1;
        b = round(3/pi*sqrt(sigma/2)*length(x));
    end
    if isfield(params,'fmax')
        fmax = params.fmax;
    else
        fmax = 0.5;
    end
    if isfield(params,'r_max')
        r_max = params.r_max;
    else
        r_max = 0;
    end
    if isfield(params,'Criteria')
        if isa(params.Criteria,'char')
            crit = {params.Criteria};
        else
           crit = params.Criteria;
        end
        crit_params = params.crit_params;
    else
        crit = {'Wang'};
        crit_params = struct('c',[4,6,8,10,12]);
    end
end

N = length(x);

[F,sF] = STFT_Gauss(x,N,sigma,fmax);

c = ridge_ext(F,0.1,0.1,10,10);
if r_max == 0
    r_max = floor(0.5*N/max(c));
end

if (max(c)<0.5*N)&&(sum(c==0)==0)
    [A1,phi1] = extract_harmonics(F,sF,c,b,b,1);

    r_opt = order_opt(x,r_max,A1,phi1,crit,crit_params);
    if r_opt>r_max
        r_opt = r_max;
    end

    [AD,phiD] = extract_harmonics(F,sF,c,b,b,r_opt);
else
    AD = nan;
    phiD = nan;
end

    TFdata.F = F;
    TFdata.sigma = sigma;
    TFdata.fmax = fmax;
    TFdata.c = c;
    TFdata.b = b;
end