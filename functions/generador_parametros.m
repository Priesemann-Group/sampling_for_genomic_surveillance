function [Par,Fun] = generador_parametros(setup)

%% Parametros

Par.k       = setup.k;
Par.gamma   = 0.1;
Par.M       = 1e6;
Par.epsilon = 5e-2;
Par.nu      = 5e-2;
Par.rho     = 1/4;
Par.xi      = setup.xi;
Par.eta     = setup.eta;
Par.lambda_s= 0.25;
Par.lambda_r= 0;
Par.R0max   = setup.R0max;
Par.ID      = generador_indices(Par.k);
Par.cutOff  = 0;
Par.cota_min  = 0.01;
Par.min_time = 1*7;
Par.Phibase = 10;
Par.tspan = setup.tspan;

%% Generador de condiciones aleatorias controladas

Par.nvars = 1;
Par.isVarThere  = zeros(Par.k,1);
Par.isVarThere(1:Par.nvars) = 1;
Par.I0 = 1000*Par.isVarThere;

%% Diseño Influx

Par.Phimax  = setup.Phimax;
Par.a       = setup.a;
Par.b       = setup.b;
Par.maxvar  = setup.maxvar;
Par.Tin     = setup.Tin;
for i = 1:Par.k
    Par.Phimax(i) = Par.Phimax(i);%/gampdf(min(Par.a(i)*Par.b(i)^2,...
%         Par.maxvar),Par.a(i),Par.b(i));
end

%% Funciones

Fun.Rt      = @(t,Par) Rt(t,Par);
Fun.Phi     = @(t,Par) Phi(t,Par);
