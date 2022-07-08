%%% main extraction statistical results into csv files %%%

clear all
close all
clc

%% Defining loop parameters

Kindice = [25 50 100];

for expno = 1:4
    for i = 1:3
        P.K = Kindice(i);
%         Matrix{expno,i} = nan(m_max,2*(k-1));
        load(strcat('Fig_1_Experiment_', num2str(expno),'_K_',num2str(P.K),'_mode_0.mat'),'Tdet_var')
        filename = strcat('Const_Fig_1_Exp_',num2str(expno),'_K_',num2str(P.K),'.csv');
        T = array2table(Tdet_var(:,2:end));
        T.Properties.VariableNames(1:4) = {'Variant_2','Variant_3','Variant_4','Variant_5'};
        writetable(T,filename)
        load(strcat('Fig_1_Experiment_', num2str(expno),'_K_',num2str(P.K),'_mode_1.mat'),'Tdet_var')
        filename = strcat('Adapt_Fig_1_Exp_',num2str(expno),'_K_',num2str(P.K),'.csv');
        T = array2table(Tdet_var(:,2:end));
        T.Properties.VariableNames(1:4) = {'Variant_2','Variant_3','Variant_4','Variant_5'};
        writetable(T,filename)
    end
end

