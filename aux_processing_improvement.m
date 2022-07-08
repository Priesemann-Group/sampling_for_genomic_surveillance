%%% quantifying improvements due to sampling protocols %%%

clear all
close all
clc

%% Importing relevant variables

load('Ground_truth.mat')
load('DefColors.mat')
load('Colores.mat')
load('Experimento_Definitivo_mode0.mat')
Mediana_const = nan(k-1,length(K_ind));
Mediana_adapt = nan(k-1,length(K_ind));

Tin_true = zeros(1,4);
for i = 1:k-1
    [~,B] = find(Nobs_com_weekly(i+1,:)>0,1,'first');
    Tin_true(i) = T(B);
end

for i = 1:k-1
    for j = 1:length(K_ind)
        Mediana_const(i,j) = quantile(Tdet_com{i}(:,j)-Tin_true(i),.5);
    end
end

load('Experimento_Definitivo_mode1.mat')

Tin_true = zeros(1,4);
for i = 1:k-1
    [~,B] = find(Nobs_com_weekly(i+1,:)>0,1,'first');
    Tin_true(i) = T(B);
end

for i = 1:k-1
    for j = 1:length(K_ind)
        Mediana_adapt(i,j) = quantile(Tdet_com{i}(:,j)-Tin_true(i),.5);
    end
end

%% quantifying the improvement

K0 = 25;
median = Mediana_adapt(4,:);
K = K_ind;
fun_delay =  @(x) x(1).*exp(-(K-K0)./x(2)) + x(3);
fun_sq = @(x) sum((fun_delay(x) - median).^2);
x = fmincon(fun_sq,[12 50 0],[],[],[],[],[2 K0 0],[30 400   5]);

h = figure('units','centimeters','position',[3,3,5,5]);

plot(K_ind,median)
hold on
plot(K_ind,fun_delay(x),'--')
ylim([0 20])



% Improvement = (Mediana_const);
Improvement = (Mediana_adapt);
% Improvement = (Mediana_const);

%% figura simple


h = figure('units','centimeters','position',[3,3,15,5]);
for i = 1:k-1
    subplot(1,4,i)
    plot(K_ind,(Mediana_const(i,:)-Mediana_adapt(i,:))./Mediana_const(i,:),'Color',Default(i+1,:),'LineWidth',2)
    hold on
    xlabel('Total seq K','FontSize',10)
    title(strcat('GS mode ',num2str(mode)),'FontSize',12)
    if i == 1
        ylabel('Improvement [%]','FontSize',10)
    end
    xlim([K_ind(1) K_ind(end)])
    set(gca,'FontSize',8,'XColor','k','YColor','k')
end

%

% print(h,strcat('GS_mode_',num2str(mode),'_kcom_',num2str(P.K)),'-dpdf')
% print(h,strcat('GS_mode_',num2str(mode),'_kcom_',num2str(P.K)),'-djpeg')

