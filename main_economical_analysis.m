%%% main econ analysis %%%

clear all
close all
clc

%% Importing variabels

load('Ground_truth.mat')
load('DefColors.mat')
load('Colores.mat')

%% Parameters for (marginal) cost and utility functions

BC = 1e2;
D0 = 15;
K0 = 25;

Klim1 = 200;
Klim2 = 300;
alpha1 = 5;
alpha2 = 2.5;

m = 5e3;
alpha = 5e2;
beta = 0.75;
Kref  = 50;
D = @(K) D0 * exp(-K./Kref);
d = @(K) D0-D(K);

C = @(K) m*K + BC;
U = @(K) alpha*d(K).^beta; 

MC = @(K) (m/(D0*Kref)) * exp(K./Kref);
MU = @(K) alpha*beta*d(K).^(beta-1);

dCdK = @(K) m*ones(size(K));
dUdK = @(K) alpha * beta * ((D0^(beta))/Kref) * (1-exp(-(K-K0)./Kref)).^(beta-1) .* exp(-(K-K0)./Kref);

K = linspace(K0,5e2,1e3);
K1 = linspace(K0,Klim1-1e-2,1e3);
K2 = linspace(K0,Klim2-1e-2,1e3);
h = figure('units','centimeters','position',[3,3,5,5]);

plot(K,dUdK(K),'LineWidth',2)
hold on
% yyaxis right
plot(K1,alpha1*K1./(Klim1-K1),'LineWidth',2)
hold on
plot(K2,alpha2*K2./(Klim2-K2),'LineWidth',2)
hold on
ylim([0 50])
% plot(K,dCdK(K))
% hold on
% plot(K,dUdK(K))

% ylim([0 100])
xlim([K0 325])
set(gca,'FontSize',8,'XColor','k','YColor','k')%,'TickLength',[0.025 0.025])
xlabel('Total seq K','FontSize',10)
ylabel('economical equivalents','FontSize',10)

%% Importing previous results

%load('Experimento_Definitivo_mode1.mat')
% load('results_simulation.mat')

Tin_true = zeros(1,4);
for i = 1:k-1
    [~,B] = find(Nobs_com_weekly(i+1,:)>0,1,'first');
    Tin_true(i) = T(B);
end

Q_com = cell(k-1,1);
Q_POE = cell(k-1,1);

for i = 1:k-1
    Q_com{i} = nan(length(q),length(K_ind));
    Q_POE{i} = nan(length(q),length(K_ind));
    for j = 1:length(K_ind)
        Q_com{i}(:,j) = quantile(Tdet_com{i}(:,j)-Tin_true(i),q);
        Q_POE{i}(:,j) = quantile(Tdet_POE{i}(:,j)-Tin_true(i),q);
    end
end

h = figure('units','centimeters','position',[3,3,18,18]);
for i = 1:k-1
    subplot(2,2,i)
    
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
    % Plotting the upper interval for the 68% CI
    patch([K_ind fliplr(K_ind)],[ysim_q(3,:) fliplr(ysim_q(4,:))],Default(i+1,:),'FaceAlpha',0.25,'LineStyle','none','HandleVisibility','off');%,[0 0.5 0],'FaceColor','interpolate')
    hold on
    % Plotting the median trend
    plot(K_ind,ysim_q(3,:),'Color',Default(i+1,:),'LineWidth',2)
    hold on
    xlabel('Total number of sequenced samples K','FontSize',10)
    title(strcat('GS mode ',num2str(mode)),'FontSize',12)
    ylabel('Delay in detection','FontSize',10)
    xlim([K_ind(1) K_ind(end)])
    ylim([0 10])
    xlim([0 250])
    set(gca,'FontSize',8,'XColor','k','YColor','k')
end

%

% print(h,strcat('GS_mode_',num2str(mode),'_kcom_',num2str(P.K)),'-dpdf')
% print(h,strcat('GS_mode_',num2str(mode),'_kcom_',num2str(P.K)),'-djpeg')

