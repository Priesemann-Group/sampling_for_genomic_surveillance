%%% main expfit %%%

clear all
close all
clc

%% Importing data and relevant variables

load('Datos_para_expfit.mat')

%% Fitting parameters

Mode = [0 1];
K0 = K_ind(1);
K = K_ind;
fun_delay =  @(x,K,K0) x(1).*exp(-(K-K0)./x(2)) + x(3);
params_expfit = cell(4,4,2);
for I = 1:2
    for expno = 1:4
        for i = 1:4
            median = Medianas{expno,i,I};
            fun_sq = @(x) sum((fun_delay(x,K,K0) - median).^2);
            params_expfit{expno,i,I} = fmincon(fun_sq,[12 50 0],[],[],[],[],[2 K0 0],[30 400 5]);
        end
    end
end
save('Params_expfit.mat','params_expfit','fun_delay')


