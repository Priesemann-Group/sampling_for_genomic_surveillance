function [k,setup] = SetupExperiment(expno,tspan)

if expno == 1
    k           = 5;
    setup.k     = k;          % Número de variantes en competencia (todo debe calzar)
    setup.xi    = 0.25*ones(setup.k,1);
    setup.eta   = 0.95*ones(setup.k,1);
    setup.R0max = [1.5 1.5 1.8 2 3]';
    setup.a     = [3 4 3 3 2];
    setup.b     = [4 5 6 5 5];
    setup.maxvar= 6*30;
    setup.Phimax= [2 3 2 3 1]*5e2;
    setup.Tin   = [-5 50 75 125 200];
    setup.tspan = tspan;
elseif expno == 2
    k           = 5;
    setup.k     = k;          % Número de variantes en competencia (todo debe calzar)
    setup.xi    = 0.25*ones(setup.k,1);
    setup.eta   = 0.95*ones(setup.k,1);
    setup.R0max = [1.6 1.8 2 2.5 3]';
    setup.a     = [3 4 3 3 2];
    setup.b     = [4 5 6 5 5];
    setup.maxvar= 6*30;
    setup.Phimax= [2 3 2 3 1]*5e2;
    setup.Tin   = [-5 30 100 150 150];
    setup.tspan = tspan;

elseif expno == 3
    k           = 5;
    setup.k     = k;          % Número de variantes en competencia (todo debe calzar)
    setup.xi    = 0.25*ones(setup.k,1);
    setup.eta   = 0.95*ones(setup.k,1);
    setup.R0max = [1.5 1.75 2 2.25 4.5]';
    setup.a     = [3 4 3 3 2];
    setup.b     = [4 5 6 5 5];
    setup.maxvar= 6*30;
    setup.Phimax= [2 3 2 3 1]*5e2;
    setup.Tin   = [-5 75 100 125 250];
    setup.tspan = tspan;
elseif expno == 4
    k           = 5;
    setup.k     = k;          % Número de variantes en competencia (todo debe calzar)
    setup.xi    = 0.25*ones(setup.k,1);
    setup.eta   = 0.95*ones(setup.k,1);
    setup.R0max = [1.5 1.75 2 2.25 3.5]';
    setup.a     = [3 4 3 3 2];
    setup.b     = [4 5 6 5 5];
    setup.maxvar= 6*30;
    setup.Phimax= [2 3 2 3 1]*5e2;
    setup.Tin   = [-5 75 100 125 250];
    setup.tspan = tspan;
else
    error('Experiment undefined')
end