function [kcom,kPOE] = adapt_kcom(T,i,V_sampli_com,V_sampli_POE,m,kcom,kPOE,Nobs_POE_weekly,params)

k = params.k;
NTPOE = sum(Nobs_POE_weekly);
if NTPOE > 0
    fi_hat_com = zeros(k,i);
    fi_hat_POE = zeros(k,i);
    
    for j = 1:k
        fi_hat_com(j,:) = V_sampli_com{j}(m,1:i);
        fi_hat_POE(j,:) = V_sampli_POE{j}(m,1:i);
    end
    
    % Calculo Lambda
    Lambda = max(fi_hat_POE(:,end-1) - fi_hat_com(:,end-1));

    % Calculo Theta
    Theta = 100*max(fi_hat_com(:,end-1)-fi_hat_com(:,end-2));
    % calculo alpha
    
    alpha = 1/(1+exp(-params.s_alpha*(Lambda-params.Lambda12)));
    zeta = 1/(1+exp(params.s_zeta*(100*Theta-params.Theta12)));
    f = (params.f_alpha*alpha-zeta);
    f_norm = f/params.f_alpha;                                                             % Posibilidad de escalar con más parámetros arbitrarios...
    kcom(m,i) = kcom(m,i-1) + floor(f_norm*params.Deltakcom_max);

    if kcom(m,i) > params.kcom_max
        kcom(m,i) = params.kcom_max;
    elseif kcom(m,i)< params.kcom_min
        kcom(m,i) = params.kcom_min;
    end
    
    kPOE(m,i) = params.K-kcom(m,i);
else
    kPOE(m,i) = NTPOE;
    kcom(m,i) = params.K-kPOE(m,i);
end