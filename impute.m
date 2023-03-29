function [x_imp] = impute(x,st,L,method,params)
%%
% Wrapper function for imputation methods. Can be used to run 
% multiple imputation methods at once. Returns a matrix with each row
% containing an imputed signal.

% Inputs:
%         x: signal with missing data.
%         st: initial indexes of missing data intervals.
%         L: lengths of missing data intervals.     
%         method: either a cell of strings with the names of the
%                 built-in methods or a function handle to a custom
%                 imputation method. The custom imputation method has 
%                 to follow the standard definition: 
%                   function [s_imp] = method(s,st,L,params)
%         params: optional parameter struct or cell of structs.
% Outputs:
%       x_imp: signals with imputed values.

if nargin<5
    params = '0';
end
if nargin<4
    x_imp = impute_tlm(x,st,L);
else
    if iscell(method)
        if isa(params,'struct')
            params = {params};
        end
        
        K = size(method,2);
        x_imp = zeros(K,length(x));
        for k=1:K
            switch lower(method{k})
                case 'tlm'
                    if nargin<5
                        x_imp(k,:) = impute_tlm(x,st,L);
                    else
                        
                        x_imp(k,:) = impute_tlm(x,st,L,params{k});
                    end
                case 'lse'
                    if nargin<5
                        x_imp(k,:) = impute_lse(x,st,L);
                    else
                        
                        x_imp(k,:) = impute_lse(x,st,L,params{k});
                    end
                case 'dmd'
                    if nargin<5
                        x_imp(k,:) = impute_dmd(x,st,L);
                    else
                        
                        x_imp(k,:) = impute_dmd(x,st,L,params{k});
                    end
                case 'gpr'
                    if nargin<5
                        x_imp(k,:) = impute_gpr(x,st,L);
                    else
                        
                        x_imp(k,:) = impute_gpr(x,st,L,params{k});
                    end
                case 'arimaf'
                    if nargin<5
                        x_imp(k,:) = impute_arimaf(x,st,L);
                    else
                        
                        x_imp(k,:) = impute_arimaf(x,st,L,params{k});
                    end
                case 'arimab'
                    if nargin<5
                        x_imp(k,:) = impute_arimab(x,st,L);
                    else
                        
                        x_imp(k,:) = impute_arimab(x,st,L,params{k});
                    end
                otherwise
                    error(['Unknown imputation method: ' method{k}])
            end
            
        end
    else
        if nargin<5
            x_imp = method(x,st,L);
        else
            x_imp = method(x,st,L,params);
        end
    end
end
