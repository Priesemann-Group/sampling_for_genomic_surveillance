function [D_com, D_POE, E_com, E_POE] =  detect_delay_and_mismatch(t,y,Par,Fun,setup,seq_rate,mode,m_max)
% valid only for two variants

k = setup.k;
ti = setup.ti;
tf = setup.tf;
P.K = seq_rate;
[Nobs_com, Nobs_POE, N_com, N_POE] = Nobs(t,y,Par,Fun);
q = 0.5;      % median 

%% Parametros C(A)RD

P.kcom_i = floor(0.6*P.K);
P.kcom_max = floor(0.95*P.K);
P.kcom_min = P.kcom_i;
P.k = k;
P.s_alpha = 5;
P.s_zeta = 1;
P.Lambda12 = .25;
P.Theta12 = 10;
P.Deltakcom_max = 10;
P.f_alpha = 1.25;


%% Accumulating and discretising weekly cases

Nobs_com_weekly = nan(k,floor((tf-ti)/7) + 1);
Nobs_POE_weekly = nan(k,floor((tf-ti)/7) + 1);
T = floor(ti/7):floor(tf/7);
for i = 2:floor((tf-ti)/7) + 1
    for j = 1:k
        idx = t/7 >= T(i-1) & t/7 <= T(i);
        Nobs_com_weekly(j,i) = floor(trapz(t(idx),Nobs_com(j,idx)));
        Nobs_POE_weekly(j,i) = floor(trapz(t(idx),Nobs_POE(j,idx)));
    end
end

%% Calculating the true time of introduction

Tin_true = zeros(1,k-1);
for i = 1:k-1
    [~,B] = find(Nobs_com_weekly(i+1,:)>0,1,'first');
    Tin_true(i) = T(B);
end

%% sorting and sampling in weekly cases


T(1) = [];
Nobs_com_weekly(:,1) = [];
Nobs_POE_weekly(:,1) = [];

kcom = P.kcom_i*ones(size(T));
kcom = repmat(kcom,m_max,1);
kPOE = P.K - kcom;
V_sampli_com = cell(k,1);
V_sampli_POE = cell(k,1);
Det_delay_com = cell(4,1);
Det_delay_POE = cell(4,1);

for i = 1:k-1
    Det_delay_com{i} = NaN(m_max,1);
    Det_delay_POE{i} = NaN(m_max,1);
end


if mode == 0
    for variant = 1:k
        V_aux_com = nan(m_max,floor((tf-ti)/7));
        V_aux_POE = nan(m_max,floor((tf-ti)/7));
        for m = 1:m_max
            rng(m)
            for i = 1:floor((tf-ti)/7)
                % community
                fi_hat_com = GS(Nobs_com_weekly(:,i),kcom(m,i),m);
                V_aux_com(m,i) = fi_hat_com(variant);
                % POE
                fi_hat_POE = GS(Nobs_POE_weekly(:,i),kPOE(m,i),m);
                V_aux_POE(m,i) = fi_hat_POE(variant);
            end
        end
        V_sampli_com{variant} = V_aux_com;
        V_sampli_POE{variant} = V_aux_POE;
    end
    
    for variant = 1:k-1
        Tdet_var_com = nan(m_max,1);
        Tdet_var_POE = nan(m_max,1);
        for m = 1:m_max
            try
                [~,B] = find(V_sampli_com{variant+1}(m,:)>0,1,'first');
                Tdet_var_com(m) = T(B);
            catch err
            end
            try
                [~,B] = find(V_sampli_POE{variant+1}(m,:)>0,1,'first');
                Tdet_var_POE(m) = T(B);
            catch err
            end
        end
        Det_delay_com{variant} = Tdet_var_com -Tin_true(variant);
        Det_delay_POE{variant} = Tdet_var_POE -Tin_true(variant);
    end
else
    for i = 1:k
        V_sampli_com{i} = nan(m_max,floor((tf-ti)/7));
        V_sampli_POE{i} = nan(m_max,floor((tf-ti)/7));
    end
    
    for i = 1:floor((tf-ti)/7)
        for m = 1:m_max
            if i>2
                [kcom,kPOE] = adapt_kcom(T,i,V_sampli_com,V_sampli_POE,m,kcom,kPOE,Nobs_POE_weekly(:,i),P);
            end
            fi_hat_com = GS(Nobs_com_weekly(:,i),kcom(m,i),m);
            fi_hat_POE = GS(Nobs_POE_weekly(:,i),kPOE(m,i),m);
            for variant = 1:k
                V_sampli_POE{variant}(m,i) = fi_hat_POE(variant);
                V_sampli_com{variant}(m,i) = fi_hat_com(variant);
            end
        end
    end
    for variant = 1:k-1
        Tdet_var_com = nan(m_max,1);
        Tdet_var_POE = nan(m_max,1);
        for m = 1:m_max
            try
                [~,B] = find(V_sampli_com{variant+1}(m,:)>0,1,'first');
                Tdet_var_com(m) = T(B);
            catch
            end
            [~,B] = find(V_sampli_POE{variant+1}(m,:)>0,1,'first');
            try
                Tdet_var_POE(m) = T(B);
            catch
            end
        end
        Det_delay_com{variant} = Tdet_var_com -Tin_true(variant);
        Det_delay_POE{variant} = Tdet_var_POE -Tin_true(variant);
    end
end


t_week = t'/7;
fi_com = N_com./sum(N_com);
fi_POE = N_POE./sum(N_POE);
fi_gt_com = zeros(k,length(T));
fi_gt_POE = zeros(k,length(T));

for i = 1:k
    fi_gt_com(i,:) = linterp(t_week,fi_com(i,:),T);
    fi_gt_POE(i,:) = linterp(t_week,fi_POE(i,:),T);
end

D_com = Det_delay_com{1};
D_POE = Det_delay_POE{1};
E_com = V_sampli_com{1}-fi_gt_com(2,:);
E_POE = V_sampli_POE{1}-fi_gt_com(2,:);

