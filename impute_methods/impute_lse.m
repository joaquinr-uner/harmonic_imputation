function [s_imp] = impute_lse(s,st,L,params)
%%
% Missing data imputation by forecasting using lease square estimation of 
% the linearized dynamical model.
% Inputs: 
%         s: signal with missing data.
%         st: initial indexes of missing data intervals.
%         L: lengths of missing data intervals.         
%         params: optional parameters struct. Includes:
%                  K: size of the dataset for forecasting
%                  M: lengths of the segments used for forecasting
% Outputs:
%         s_imp: signal with imputed values.
if nargin<4
    K = 0;
    M = 0;
    %M = round(1*L);
    %K = round(1.5*M);
else
    if isfield(params,'M')
        M = params.M;
        if length(M)==1
            M = M*ones(1,length(L));
        end
    else
        M = zeros(1,length(L));
    end
    if isfield(params,'K')
        K = params.K;
        if length(K)==1
            K = K*ones(1,length(L));
        end
    else
        K = round(2.5*M);
    end
end

s = s(:);

N = length(s);
Ni = length(st);

s_imp = s;
intp = 0;
for qi=1:Ni
    inti = st(qi):st(qi)+L(qi)-1;
    sp = s_imp(intp(end)+1:inti(1)-1);
    Np = length(sp);
    Ki = K(qi);
    Mi = M(qi);
    if Ki==0
        [~,fh] = compute_sigma(sp,1);
        Th = length(sp)/fh;
        %Mi = floor(0.9*Th);
        %Ki = ceil(3*Th);
        Mi = ceil(3*Th);
        Ki = floor(2.5*Mi);

    end

    if Ki+Mi>Np
        Ki = floor(0.7*Np);
        Mi = floor(0.3*Np)-1;
    end
    %sp = s_imp(1:inti(1)-1);

    X = zeros(Mi,Ki) ;
    for k = 1: Ki
        X(:,k) = sp(end-Ki-Mi+k:end-Ki+k-1);
    end
    Y = [X(:,2:end) sp(end-Mi+1:end)] ;

    A = (Y*X') / (X*X') ; % least square estimation


    Z = zeros(Mi,L(qi)) ;
    Z(:,1) = A*Y(:,end) ;
    for k = 2:L(qi)
        Z(:,k) = A*Z(:,k-1) ;
    end
    xext = Z(end,:)' ;

    s_imp(inti) = xext;
    intp = inti;
end
