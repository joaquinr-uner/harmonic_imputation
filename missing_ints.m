function [st,LEN,Ni] = missing_ints(s, d, t)
%%
% Find missing data intervals in signal s, with threshold t and minimum
% interval length d.
% Inputs: 
%         s: signal with missing data.
%         t: threshold.
%         d: minimun missing data interval.
% Outputs:
%         st: initial indexes of missing data intervals.
%         LEN: lengths of missing data intervals.
%         Ni: number of distinct missing data intervals.

    tmp = abs(s)<=t;
    [M, st] = regexp(sprintf('%i', tmp'), '1+', 'match') ;
    
    p = 0;
    for i = 1:size(M,2)
        if length(M{i})>d
            p = p + 1;
            LEN(p) = length(M{i});
            aux(p) = st(i);
        end
    end
    st = aux;
    Ni = length(LEN);
end