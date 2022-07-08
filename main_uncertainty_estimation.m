%%% main uncertainty estimator %%%

clear all
close all
clc

%% Defining loop and loop parameters

Mode = [0 1];
funError = @(M,v) abs(v-M);%./(v);
sp_id = [(1:2:10)' 1+(1:2:10)'];

for expno = 1:4
    for Kexpeno = [25 50 100]
        h = figure('units','centimeters','position',[15,3,6,9]);
        for i_mode = 1:2
            mode = Mode(i_mode);
            P.K = Kexpeno;
            filename = strcat('Fig_1_Experiment_', num2str(expno),'_K_',num2str(P.K),'_mode_',num2str(mode),'.mat');
            load(filename)
            
            Q = cell(k,1);
            for i = 1:k
                Q{i} = nan(length(q),floor((tf-ti)/7));
                for j = 1:floor((tf-ti)/7)
                    Q{i}(:,j) = quantile(funError(V_sampli_com{i}(:,j),fi(i,j)),q);
                end
            end
            
            Taux = T;
            for i = 1:k
                idx = fi(i,:)>0.01;
                T = Taux(idx);
                % Plotting the lower interval for the 95% CI
                subplot(5,2,sp_id(i,i_mode))
                ysim_q = Q{i}(:,idx);%./fi(i,:);
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
                % Plotting the median trend
                plot(T,ysim_q(3,:),'Color',Default(i,:),'LineWidth',2*fact_curva)
                hold on
                ylim([0 0.25])
                xlim([Taux(1) Taux(end)])
            end
        end
        set(gca,'FontSize',8,'XColor','k','YColor','k','TickLength',[0.025 0.025])
        print(h,strcat('Sfig_var_Exp_no_',num2str(expno),'_K_',num2str(P.K)),'-dpdf')
        close all
        clear h
    end
end
