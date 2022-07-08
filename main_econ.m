%%% main figure Kobs b, worldwide %%%

clear all
close all
clc

%% Importing data

load('Obs_Seq_Data.mat')

aux = unique(Country);
% aux = aux(1:end-1);
inHabs = [45.38 212.6 19.12 5.83 83.24 5 59.55 17.44 32.97 37.95 59.31 47.35 67.22 329.5];
max_seq = zeros(size(aux));
for i = 1:length(aux)
    idx = Country == aux(i);
    freq_seq(idx) = freq_seq(idx)/inHabs(i);
    max_seq(i) = max(freq_seq(idx));
end
ub = [50 300 500 1200];
lb = [0 50 300 500];

%% plotting

h = figure('units','centimeters','position',[3,3,15,8]);
a = boxplot(freq_seq,Country,'PlotStyle','compact','Color','k');
set(gca, 'YScale','Log')
ylim([0.1 5e3])

% for i = 1:length(ub)
%     subplot(1,4,i)
%     idx = max_seq >= lb(i) & max_seq <= ub(i);
%     idx_1 = ismember(Country,aux(idx));
%     a = boxplot(freq_seq(idx_1),Country(idx_1),'PlotStyle','compact','Color','k')
%     set(a,{'linew'},{1.5})
set(gca,'FontSize',8,'XColor','k','YColor','k','TickLength',[0.025 0.025])
%     ylim([0 ub(i)])
%     hold on
% end