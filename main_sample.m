% function h =  main_sample(expno,mode)
%mode 0 = constant; mode 1 = adaptive
for expno = 1
    for mode = 0:1
        for Kindice = [25 50 100]
            P.K = Kindice;
            %% Importing variables
            
            load(strcat('Ground_truth_Experiment_', num2str(expno),'.mat'))
            load('DefColors.mat')
            load('Colores.mat')
            
            m_max = 100;               % number of random seeds to create the fi trends;
            q = [.025 .16 0.5 0.84 .975];
            labels_resorted = {'variant 2','variant 3','variant 4','variant 5'};
            
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
            
            % fig = figure('units','cm','position',[0,0,9,9]);
            
            h = figure('units','centimeters','position',[3,3,6,9]);
            subplot(2,1,1)
            plot(T,Nobs_com_weekly,'LineWidth',2)
            hold on
            N_T = sum(Nobs_com_weekly);
            plot(T,N_T,'k-','LineWidth',2)
            xlim([0 T(end)])
            ylabel('Weekly new cases','FontSize',10)
            set(gca,'FontSize',8,'XColor','k','YColor','k','TickLength',[0.025 0.025])
            subplot(2,1,2)
            tau = 4;
            NT_cont = sum(Nobs_com);
            NT_cont_daily = interp1(t,NT_cont,ti:tf);
            
            
            Rt_hat = NT_cont_daily./[NaN(1,tau) NT_cont_daily(1:end-tau)];
            plot((ti:tf)/7,Rt_hat,'k-','LineWidth',2)
            hold on
            plot([0 T(end)],[1 1],'r--','LineWidth',1)
            ylim([0.8 1.2])
            xlim([0 T(end)])
            ylabel('obs. Rep number','FontSize',10)
            xlabel('weeks','FontSize',10)
            set(gca,'FontSize',8,'XColor','k','YColor','k','TickLength',[0.025 0.025])
            
%             print(h,strcat('Experiment_no_',num2str(expno)),'-dpdf')
            %print(h,strcat('GroundTruth_',num2str(mode),'_kcom_',num2str(P.K)),'-djpeg')
            
            idx = sum(Nobs_POE_weekly) == 0;        % indices where no cases were detected at POES
            
            %% sorting and sampling from weekly cases
            
            T(1) = [];
            Nobs_com_weekly(:,1) = [];
            Nobs_POE_weekly(:,1) = [];
            
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
            end
            
            Q = cell(k,1);
            for i = 1:k
                Q{i} = nan(length(q),floor((tf-ti)/7));
                for j = 1:floor((tf-ti)/7)
                    Q{i}(:,j) = quantile(V_sampli_com{i}(:,j),q);
                end
            end

            
            h = figure('units','centimeters','position',[3,3,6,9]);
            subplot(3,1,2)
            for i = 1:k
                % Plotting the lower interval for the 95% CI
                ysim_q = Q{i};
                patch([T fliplr(T)],[ysim_q(1,:) fliplr(ysim_q(2,:))],Default(i,:),'FaceAlpha',0.15,'LineStyle','none','HandleVisibility','off');%,[0 0.5 0],'FaceColor','interpolate')
                hold on
                % Plotting the upper interval for the 95% CI
                patch([T fliplr(T)],[ysim_q(4,:) fliplr(ysim_q(5,:))],Default(i,:),'FaceAlpha',0.15,'LineStyle','none','HandleVisibility','off');%,[0 0.5 0],'FaceColor','interpolate')
                hold on
                % Plotting the lower interval for the 68% CI
                patch([T fliplr(T)],[ysim_q(2,:) fliplr(ysim_q(3,:))],Default(i,:),'FaceAlpha',0.25,'LineStyle','none','HandleVisibility','off');%,[0 0.5 0],'FaceColor','interpolate')
                hold on
                % Plotting the upper interval for the 68% CI
                patch([T fliplr(T)],[ysim_q(3,:) fliplr(ysim_q(4,:))],Default(i,:),'FaceAlpha',0.25,'LineStyle','none','HandleVisibility','off');%,[0 0.5 0],'FaceColor','interpolate')
                hold on
                % Plotting the median trends
                plot(T,ysim_q(3,:),'Color',Default(i,:),'LineWidth',2)
                hold on
            end
            
            xlabel('weeks','FontSize',10)
            title(strcat('GS_mode_',num2str(mode),', kcom =',num2str(P.K)),'FontSize',12)
            ylabel('variant share','FontSize',10)
            xlim([0 T(end)])
            set(gca,'FontSize',8,'XColor','k','YColor','k','TickLength',[0.025 0.025])
            
            subplot(3,1,1)
            plot(t/7,N_com./sum(N_com),'LineWidth',2)
            xlim([0 T(end)])
            title('Ground truth','FontSize',12)
            ylabel('variant share','FontSize',10)
            xlabel('weeks','FontSize',10)
            set(gca,'FontSize',8,'XColor','k','YColor','k','TickLength',[0.025 0.025])
            
            
            subplot(3,1,3)
            plot(T,kcom(1,:),'LineWidth',2)%(1,:))
            hold on
            plot(T,kPOE(1,:),'LineWidth',2)%(1,:))
            xlim([0 T(end)])
            ylim([0 P.K+1])
            ylabel('Sequenced samples','FontSize',10)
            legend({'kcom','kPOE'},'FontSize',8)
            set(gca,'FontSize',8,'XColor','k','YColor','k','TickLength',[0.025 0.025])
            
            
            %print(h,strcat('GS_mode_',num2str(mode),'_kcom_',num2str(P.K)),'-dpdf')
            %print(h,strcat('GS_mode_',num2str(mode),'_kcom_',num2str(P.K)),'-djpeg')
            
            
            %% Checking how well we estimated Ti
            
            Tin = setup.Tin(2:end)/7;
            Tdet_var = nan(m_max,k);
            h = figure('units','centimeters','position',[3,3,6,9]);
            for i = 2:k
                for m = 1:m_max
                    [~,B] = find(V_sampli_com{i}(m,:)>0,1,'first');
                    Tdet_var(m,i) = T(B);
                end
            end
            for i = 1:length(Tin)
                plot([Tin(i) Tin(i)],[0 5],'Color',Default(i+1,:),'LineWidth',2)
                hold on
            end
            a = boxplot(Tdet_var(:,2:end),labels_resorted,'PlotStyle','compact','LabelOrientation','horizontal','Orientation','horizontal','Color',Default(2:5,:));
            set(a,{'linew'},{1.5})
            set(gca,'FontSize',8,'XColor','k','YColor','k','TickLength',[0.025 0.025])
            
            xlim([0 T(end)])
            %print(h,strcat('Tdect_mode_',num2str(mode),'_K_',num2str(P.K),'_kcomi_',num2str(P.kcom_i)),'-dpdf')
            %print(h,strcat('Tdect_mode_',num2str(mode),'_K_',num2str(P.K),'_kcomi_',num2str(P.kcom_i)),'-djpeg')
            
            %% absolute error in variant share, per variant, over time
            
            t_week = t'/7;
            % idx = t_week >= T(1) & t_week <= T(end);
            % t_week = t_week(idx);
            % fi_com = N_com(:,idx)./sum(N_com(:,idx));
            fi_com = N_com./sum(N_com);
            fi = zeros(k,length(T));
            for i = 1:k
                fi(i,:) = linterp(t_week,fi_com(i,:),T);
            end
            
            
            h = figure('units','centimeters','position',[3,3,6,9]);
            
            for i = 1:k
                % Plotting the lower interval for the 95% CI
                subplot(5,1,i)
                ysim_q = (Q{i}-fi(i,:));%./fi(i,:);
                patch([T fliplr(T)],[ysim_q(1,:) fliplr(ysim_q(2,:))],Default(i,:),'FaceAlpha',0.15,'LineStyle','none','HandleVisibility','off');%,[0 0.5 0],'FaceColor','interpolate')
                hold on
                % Plotting the upper interval for the 95% CI
                patch([T fliplr(T)],[ysim_q(4,:) fliplr(ysim_q(5,:))],Default(i,:),'FaceAlpha',0.15,'LineStyle','none','HandleVisibility','off');%,[0 0.5 0],'FaceColor','interpolate')
                hold on
                % Plotting the lower interval for the 68% CI
                patch([T fliplr(T)],[ysim_q(2,:) fliplr(ysim_q(3,:))],Default(i,:),'FaceAlpha',0.25,'LineStyle','none','HandleVisibility','off');%,[0 0.5 0],'FaceColor','interpolate')
                hold on
                % Plotting the lower interval for the 68% CI
                patch([T fliplr(T)],[ysim_q(3,:) fliplr(ysim_q(4,:))],Default(i,:),'FaceAlpha',0.25,'LineStyle','none','HandleVisibility','off');%,[0 0.5 0],'FaceColor','interpolate')
                hold on
                % Plotting the median trend
                plot(T,ysim_q(3,:),'Color',Default(i,:),'LineWidth',2*fact_curva)
                hold on
                ylim([-0.3,0.3])
                xlim([T(1) T(end)])
            end
            set(gca,'FontSize',8,'XColor','k','YColor','k','TickLength',[0.025 0.025])
            close all
            clear h
%             save(strcat('Fig_1_Experiment_', num2str(expno),'_K_',num2str(P.K),'_mode_',num2str(mode),'.mat'))
        end
    end
end
%print(h,strcat('Varianza_mode_',num2str(mode),'_K_',num2str(P.K),'_kcomi_',num2str(P.kcom_i)),'-dpdf')
%print(h,strcat('Varianza_mode_',num2str(mode),'_K_',num2str(P.K),'_kcomi_',num2str(P.kcom_i)),'-djpeg')
