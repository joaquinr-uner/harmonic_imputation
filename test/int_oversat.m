function [s_int, trend_int] = int_oversat(s_imp,st,L,intp)
%%
% Missing data imputation refinement by harmonic interpolation.
% Inputs: 
%         A: harmonic amplitudes of the imputed signal.
%         phi: harmonic phases of the imputed signal.
%         st: initial indexes of missing data intervals.
%         L: lengths of missing data intervals.
%         intp: interpolation type. spline or pchip.
%         s_imp: signal with imputed values
% Outputs:
%         s_int: signal with improved missing data imputation.
%         A_int: interpolated harmonic amplitudes.
%         phi_int: interpolated harmonic phases.

if nargin< 7
    trend = zeros(size(s_imp));
end

s_imp = s_imp(:);
s_int = s_imp;
trend_int = trend;

N = length(s_imp);
Ni = length(L);


ed = st + L - 1;
for i=1:Ni
    Li = L(i);
    sti = st(i);
    edi = sti + Li - 1;
    
    if i==1
        lint = 1:st(i)-1;
    else
        lint = ed(i-1)+1:st(i)-1;
    end

    if i==Ni
        rint = ed(i)+1:N;
    else
        rint = ed(i)+1:st(i+1)-1;
    end
    xq = [lint rint];
    int = interp1(xq,s_imp(xq),sti:edi,intp);

    s_int(sti:edi) = int;

end



