function [Err,Err_i] = compute_errors(t,x,st,L,metric)
if nargin<5
    metric = 'mae';
else
    if ~iscell(metric)
        metric = {metric};
    end
end
t = t(:);
x = x(:);

K = length(L);
M = size(metric,2);
N = length(t);

Err_i = zeros(K,size(metric,1));
Err = zeros(1,size(metric,1));
for m=1:M
    metr = lower(metric{m});

    switch metr
        case 'mae'
            mae = zeros(1,K);
            for k=1:K
                Lk = L(k);
                stk = st(k);
                edk = stk + Lk - 1;
                intk = stk:edk;
                mae(k) = norm(t(intk)-x(intk),1)/Lk;
            end
            Err_i(:,m) = mae;
            Err(m) = sum(abs(t-x))/sum(L);
        case 'mse'
            mse = zeros(1,K);
            for k=1:K
                Lk = L(k);
                stk = st(k);
                edk = stk + Lk - 1;
                intk = stk:edk;
                mse(k) = norm(t(intk)-x(intk))^2/Lk;
            end
            Err_i(:,m) = mse;
            Err(m) = sum((t-x).^2)/sum(L);
        case 'rmse'
            rmse = zeros(1,K);
            for k=1:K
                Lk = L(k);
                stk = st(k);
                edk = stk + Lk - 1;
                intk = stk:edk;
                rmse(k) = norm(t(intk)-x(intk))/Lk;
            end
            Err_i(:,m) = rmse;
            Err(m) = sqrt(sum((t-x).^2)/sum(L));
        case 'sim'
            sim = zeros(1,K);
            for k=1:K
                Lk = L(k);
                stk = st(k);
                edk = stk + Lk - 1;
                intk = stk:edk;
                rk = max(t)-min(t);
                sim(k) = sum(rk./(rk+abs(x(intk)-t(intk))))/Lk;

            end
            Err_i(:,m) = sim;
            Err(m) = sum(rk./(rk+abs(x-t)))/sum(L);

        otherwise
            error(['Unknown Performance Metric "' metr '"'])
    end
end