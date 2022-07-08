%%% plotter figura experimentos %%%

clear all
close all
clc

h = figure('units','centimeters','position',[3,3,16,4.5]);

for expno = 1:4
    load(strcat('Fig_1_Experiment_', num2str(expno),'_K_25_mode_1.mat'))
    subplot(1,4,expno)
    plot(T,Nobs_com_weekly,'LineWidth',2)
    hold on
    N_T = sum(Nobs_com_weekly);
    plot(T,N_T,'k--','LineWidth',2)
    xlim([0 T(end)])
    if expno == 1
        ylabel('weekly new cases','FontSize',10)
    end
    ylim([0 6e3])
    xlabel('calendar week','FontSize',10)
    title(strcat('Experiment no ',num2str(expno)),'FontSize',12)
    set(gca,'FontSize',8,'XColor','k','YColor','k','TickLength',[0.025 0.025])
    hold on
end

print(h,'Figure_experiments','-dpdf')