%%% direct problem solver: script to create wave patterns %%%

clear all
close all
clc

%% Defininig solver variables

expno       = 1;
ti          = -14;
tf          = 500;
tspan       = [ti tf];

%% Loading experiment-specific parameters

[k,setup]   = SetupExperiment(expno,tspan);
[Par,Fun]   = generador_parametros(setup);


%% Solver

AbsTol = 1e-12;
RelTol = 1e-9;
options= odeset('AbsTol',AbsTol,'RelTol',RelTol,'NonNegative',1:(5*k+2),'Events',@(t,y) Extinctions(t,y,Par));

odefun = @(t,y) dydt(t,y,Par,Fun);
generador_cond_inicial(Par)

[t,y] = ...
    ode45(@(t,y) odefun(t,y), tspan, generador_cond_inicial(Par), options);

[Nobs_com, Nobs_POE, N_com, N_POE] = Nobs(t,y,Par,Fun);

%% preliminary visualization
plot(t,Nobs_POE)
xlim([0 tf])
figure
plot(t,Nobs_com)
hold on
plot(t,sum(Nobs_com),'k--')
xlim([0 tf])

% save(strcat('Ground_truth_Experiment_', num2str(expno),'.mat'))

