%%% main K-swipe %%%

clear all
close all
clc

%% Definition of general variables

load('DefColors.mat')
load('Colores.mat')
mode = 0;
m_max = 100;                % number of random seeds to create the fi trends;
K_ind = 10:5:500;           % number of total samples considered throughout
q = [.025 .16 0.5 0.84 .975];
labels_resorted = {'variant 2','variant 3','variant 4','variant 5'};

for expno = 2:4
    for mode = [1 0] 
        filename = strcat('Fig_2_Exp_',num2str(expno),'_mode_',num2str(mode));        
        load(strcat('Ground_truth_Experiment_', num2str(expno),'.mat'))
        
        
        %% Empty matrices
        
        Tdet_com = cell(4,1);
        Tdet_POE = cell(4,1);
        for i = 1:k-1
            Tdet_com{i} = NaN(m_max,length(K_ind));
            Tdet_POE{i} = NaN(m_max,length(K_ind));
        end
        
        %% accumulating and discretising weekly cases
        
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
        
        T(1) = [];
        Nobs_com_weekly(:,1) = [];
        Nobs_POE_weekly(:,1) = [];
        
        %% C(A)RD parameters
        
        for ind = 1:length(K_ind)
            
            P.K = K_ind(ind);
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
            
            %% sorting and sampling in weekly cases
            
            kcom = P.kcom_i*ones(size(T));
            kcom = repmat(kcom,m_max,1);
            kPOE = P.K - kcom;
            V_sampli_com = cell(k,1);
            V_sampli_POE = cell(k,1);
            
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
                    Tdet_com{variant}(:,ind) = Tdet_var_com;
                    Tdet_POE{variant}(:,ind) = Tdet_var_POE;
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
                        % POE
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
                    Tdet_com{variant}(:,ind) = Tdet_var_com;
                    Tdet_POE{variant}(:,ind) = Tdet_var_POE;
                end
            end
        end
        
%         save(strcat(filename,'.mat'))
    end
end

