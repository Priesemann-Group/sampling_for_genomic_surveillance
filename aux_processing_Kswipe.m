%%% main sampleo y otros tratamientos %%%

clear all
close all
clc

%% Importing relevant data

load('Ground_truth.mat')
load('DefColors.mat')
load('Colores.mat')
load('Experimento_Definitivo_mode1.mat')
% load('results_simulation.mat')

Tin_true = zeros(1,4);
for i = 1:k-1
    [~,B] = find(Nobs_com_weekly(i+1,:)>0,1,'first');
    Tin_true(i) = T(B);
end

Q_com = cell(k-1,1);
Q_POE = cell(k-1,1);
Q_aux = cell(k-1,1);

for i = 1:k-1
    Q_com{i} = nan(length(q),length(K_ind));
    Q_POE{i} = nan(length(q),length(K_ind));
    Q_aux{i} = nan(length(Tdet_com{i}(:,1)),length(K_ind));
    for j = 1:length(K_ind)
        Q_com{i}(:,j) = quantile(Tdet_com{i}(:,j)-Tin_true(i),q);
        Q_POE{i}(:,j) = quantile(Tdet_POE{i}(:,j)-Tin_true(i),q);
        Q_aux{i}(:,j) = Tdet_com{i}(:,j)-Tin_true(i);
    end
end


h = figure('units','centimeters','position',[3,3,15,5]);
for i = 1:k-1
    subplot(1,4,i)
    
    % Plotting the lower interval for the 95% CI
    ysim_q = Q_com{i};
    patch([K_ind fliplr(K_ind)],[ysim_q(1,:) fliplr(ysim_q(2,:))],Default(i+1,:),'FaceAlpha',0.15,'LineStyle','none','HandleVisibility','off');%,[0 0.5 0],'FaceColor','interpolate')
    hold on
    % Plotting the upper interval for the 95% CI
    patch([K_ind fliplr(K_ind)],[ysim_q(4,:) fliplr(ysim_q(5,:))],Default(i+1,:),'FaceAlpha',0.15,'LineStyle','none','HandleVisibility','off');%,[0 0.5 0],'FaceColor','interpolate')
    hold on
    % Plotting the lower interval for the 68% CI
    patch([K_ind fliplr(K_ind)],[ysim_q(2,:) fliplr(ysim_q(3,:))],Default(i+1,:),'FaceAlpha',0.25,'LineStyle','none','HandleVisibility','off');%,[0 0.5 0],'FaceColor','interpolate')
    hold on
    % Plotting the lower interval for the 68% CI
    patch([K_ind fliplr(K_ind)],[ysim_q(3,:) fliplr(ysim_q(4,:))],Default(i+1,:),'FaceAlpha',0.25,'LineStyle','none','HandleVisibility','off');%,[0 0.5 0],'FaceColor','interpolate')
    hold on
    % Plotting the median trends
    plot(K_ind,ysim_q(3,:),'Color',Default(i+1,:),'LineWidth',2)
    hold on
    xlabel('Total seq K','FontSize',10)
    title(strcat('GS mode ',num2str(mode)),'FontSize',12)
    if i == 1
        ylabel('Delay in detection [weeks]','FontSize',10)
    end
    xlim([K_ind(1) K_ind(end)])
    ylim([0 10])
%     xlim([0 250])
    set(gca,'FontSize',8,'XColor','k','YColor','k')
end

%

% print(h,strcat('GS_mode_',num2str(mode),'_kcom_',num2str(P.K)),'-dpdf')
% print(h,strcat('GS_mode_',num2str(mode),'_kcom_',num2str(P.K)),'-djpeg')

