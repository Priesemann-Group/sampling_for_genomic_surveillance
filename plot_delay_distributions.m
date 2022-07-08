%%% main routine for plotting delay distributions %%%

clear all
close all
clc

%% Importing plotting variables

load('DefColors.mat')
load('Colores.mat')

%% Loop for importing results

mode = 0;
m_max = 100;                % number of random seeds to create the fi trends;
K_ind = 10:5:500;           % number of total samples considered throughout
q = [.025 .16 0.5 0.84 .975];
labels_resorted = {'variant 2','variant 3','variant 4','variant 5'};

Q_com = cell(4,4,2); % experiments x variants x sampling mode
Q_POE = cell(4,4,2);
Q_aux = cell(4,4,2);
Mode = [0 1];
for I = 1:2
    for expno = 1:4
        filename = strcat('Fig_2_Exp_',num2str(expno),'_mode_',num2str(Mode(I)));
        load(strcat(filename,'.mat'))
        Tin_true = zeros(1,4);
        for i = 1:k-1
            [~,B] = find(Nobs_com_weekly(i+1,:)>0,1,'first');
            Tin_true(i) = T(B);
        end
        
        for i = 1:k-1
            Q_com{expno,i,I} = nan(length(q),length(K_ind));
            Q_POE{expno,i,I} = nan(length(q),length(K_ind));
            Q_aux{expno,i,I} = nan(length(Tdet_com{i}(:,1)),length(K_ind));
            for j = 1:length(K_ind)
                Q_com{expno,i,I}(:,j) = quantile(Tdet_com{i}(:,j)-Tin_true(i),q);
                Q_POE{expno,i,I}(:,j) = quantile(Tdet_POE{i}(:,j)-Tin_true(i),q);
                Q_aux{expno,i,I}(:,j) = Tdet_com{i}(:,j)-Tin_true(i);
            end
        end
    end
end


Medianas = cell(4,4,4);
Promedios = cell(4,4,4);
load('Params_expfit.mat')
for I = 1:2
    h = figure('units','centimeters','position',[3,3,18,16]);
    for expno = 1:4
        for i = 1:k-1
            sp_id = i + (k-1)*(expno-1);
            subplot(4,k-1,sp_id)
            
            % Plotting the lower interval for the 95% CI
            ysim_q = Q_com{expno,i,I};
            patch([K_ind fliplr(K_ind)],[ysim_q(1,:) fliplr(ysim_q(2,:))],Default(i+1,:),'FaceAlpha',0.15,'LineStyle','none','HandleVisibility','off');%,[0 0.5 0],'FaceColor','interpolate')
            hold on
            % Plotting the upper interval for the 95% CI
            patch([K_ind fliplr(K_ind)],[ysim_q(4,:) fliplr(ysim_q(5,:))],Default(i+1,:),'FaceAlpha',0.15,'LineStyle','none','HandleVisibility','off');%,[0 0.5 0],'FaceColor','interpolate')
            hold on
            % Plotting the lower interval for the 68% CI
            patch([K_ind fliplr(K_ind)],[ysim_q(2,:) fliplr(ysim_q(3,:))],Default(i+1,:),'FaceAlpha',0.25,'LineStyle','none','HandleVisibility','off');%,[0 0.5 0],'FaceColor','interpolate')
            hold on
            % Plotting the upper interval for the 68% CI
            patch([K_ind fliplr(K_ind)],[ysim_q(3,:) fliplr(ysim_q(4,:))],Default(i+1,:),'FaceAlpha',0.25,'LineStyle','none','HandleVisibility','off');%,[0 0.5 0],'FaceColor','interpolate')
            hold on
            % Plotting the median trend
            Medianas{expno,i,I} = ysim_q(3,:);
            Promedios{expno,i,I} = mean(Q_aux{expno,i,I});
            plot(K_ind,ysim_q(3,:),'Color',Default(i+1,:),'LineWidth',2)
            hold on
            % Plotting the exponential fit
            plot(K_ind,fun_delay(params_expfit{expno,i,I},K_ind,K_ind(1)),'--','Color',Default(i+1,:),'LineWidth',1.5)
            
            hold on
            xlim([0 K_ind(end)])
            ylim([0 10])
            set(gca,'FontSize',8,'XColor','k','YColor','k')
        end
        print(h,strcat('Figure_3_mode_',num2str(Mode(I))),'-dpdf')
    end
end


h = figure('units','centimeters','position',[3,3,18,4]);
for expno = 1:4
    for i = 1:k-1
        subplot(1,k-1,i)
        e = (Promedios{expno,i,1}-Promedios{expno,i,2});
        col = ((expno-1)/6) * [1 1 1] +(1-(expno-1)/6) * Default(i+1,:);
        plot(K_ind,e,'Color',col,'LineWidth',2)
        hold on
        xlim([0 K_ind(end)])
        ylim([0 5])
        set(gca,'FontSize',8,'XColor','k','YColor','k')%,'TickLength',[0.025 0.025])
    end
    print(h,'Figure_4','-dpdf')
end

% save('Datos_para_expfit.mat','Medianas','K_ind')
%

% print(h,strcat('GS_mode_',num2str(mode),'_kcom_',num2str(P.K)),'-dpdf')
% print(h,strcat('GS_mode_',num2str(mode),'_kcom_',num2str(P.K)),'-djpeg')

