%%% main parameter swipe: systematically exploring varaint detection delay %%%

clear all
close all
clc

%% Loop variables

K_repeat = 10;%[25 50 100];
Mode = [0 1];
m_max = 250;
%% setting loops

for var_i = 1:length(K_repeat)
    for var_j = 1:2
        
        %% Exp cvariables
        
        mode = Mode(var_j);
        Kexp = K_repeat(var_i);
        filename = strcat('syst_swipe_mode',num2str(mode),'K',num2str(Kexp),'_250.mat');
        
        %% Defininig solver variables
        
        expno       = 1;
        setup.ti    = -14;
        setup.tf    = 500;
        tspan       = [setup.ti setup.tf];
        
        %% Defining experiment
        
        k           = 2;
        setup.k     = k;          % Número de variantes en competencia (todo debe calzar)
        setup.xi    = 0.25*ones(setup.k,1);
        setup.eta   = 0.95*ones(setup.k,1);
        setup.tspan = tspan;
        setup.a     = [3 4];
        setup.b     = [4 5];
        setup.maxvar= 6*30;
        setup.Phimax= [2 3]*5e2;
        AbsTol = 1e-12;
        RelTol = 1e-9;
        options= odeset('AbsTol',AbsTol,'RelTol',RelTol,'NonNegative',1:(5*k+2));
        
        %% calculating the reference time
        
        setup.R0max = [1.5 1.5]';
        setup.Tin   = [-5 450];
        [Par,Fun]   = generador_parametros(setup);
        odefun = @(t,y) dydt(t,y,Par,Fun);
        generador_cond_inicial(Par)
        [t,y] = ...
            ode45(@(t,y) odefun(t,y), tspan, generador_cond_inicial(Par), options);
        [~, ~, N_com, N_POE] = Nobs(t,y,Par,Fun);
        [~,ind] = max(sum(N_com));
        Tmax = 1+floor(t(ind)/7);
        
        %% we have to make a swipe
        
        R0v1 = 1.5;
        n_R0 = 21;
        R0max_swipe = R0v1*linspace(1,6,n_R0);
        R0max_swipe_ref = linspace(1,6,n_R0);
        n_Tin = 9;
        Tin_swipe = 7*(Tmax + linspace(-20,20,n_Tin));
        Tin_swipe_ref = linspace(-20,20,n_Tin);
        D_com = cell(n_R0,n_Tin);
        D_POE = cell(n_R0,n_Tin);
        E_com = cell(n_R0,n_Tin);
        E_POE = cell(n_R0,n_Tin);
        
        for I = 1:n_R0
            for J = 1:n_Tin
                setup.R0max = [R0v1 R0max_swipe(I)]';
                setup.Tin   = [-5 Tin_swipe(J)];
                [Par,Fun]   = generador_parametros(setup);
                odefun = @(t,y) dydt(t,y,Par,Fun);
                generador_cond_inicial(Par)
                [t,y] = ...
                    ode45(@(t,y) odefun(t,y), tspan, generador_cond_inicial(Par), options);
                [D_com{I,J}, D_POE{I,J}, E_com{I,J}, E_POE{I,J}] =  detect_delay_and_mismatch(t,y,Par,Fun,setup,Kexp,mode,m_max);
            end
        end
        
        
        %% saving results
        
        save(filename)
        
    end
end
%% plotting

% %contour(D_com)
% aux_D = D_com';
% aux_D1  = flip(aux_D,1);
% heatmap(round(R0max_swipe/1.5,2),round((flip(Tin_swipe)-Tmax)/7,1),aux_D1)
% xlabel('Rel. rep number second variant')
% ylabel('Int. time second variant (relative to peak) [weeks]')
% caxis([1, 10]);


