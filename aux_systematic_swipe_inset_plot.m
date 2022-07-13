%%% main parameter swipe: systematically exploring varaint detection delay %%%

clear all
close all
clc

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

Jind = [2 5 8];
Iind = [2 4 8];
ylimsup = [2500 3000 12500];
for m = 1:3
    I = Iind(m);
    J = Jind(m);
    setup.R0max = [R0v1 R0max_swipe(I)]';
    setup.Tin   = [-5 Tin_swipe(J)];
    [Par,Fun]   = generador_parametros(setup);
    odefun = @(t,y) dydt(t,y,Par,Fun);
    generador_cond_inicial(Par)
    [t,y] = ...
        ode45(@(t,y) odefun(t,y), tspan, generador_cond_inicial(Par), options);
    [Nobs_com, Nobs_POE, N_com, N_POE] = Nobs(t,y,Par,Fun);
    %% Accumulating and discretising weekly cases
    
    Nobs_com_weekly = nan(k,floor((setup.tf-setup.ti)/7) + 1);
    Nobs_POE_weekly = nan(k,floor((setup.tf-setup.ti)/7) + 1);
    T = floor(setup.ti/7):floor(setup.tf/7);
    for i = 2:floor((setup.tf-setup.ti)/7) + 1
        for j = 1:k
            idx = t/7 >= T(i-1) & t/7 <= T(i);
            Nobs_com_weekly(j,i) = floor(trapz(t(idx),Nobs_com(j,idx)));
            Nobs_POE_weekly(j,i) = floor(trapz(t(idx),Nobs_POE(j,idx)));
        end
    end
    
    
    h = figure('units','centimeters','position',[3*m,3,3,3]);
    plot(T,Nobs_com_weekly,'LineWidth',2)
    hold on
    N_T = sum(Nobs_com_weekly);
    plot(T,N_T,'k--','LineWidth',2)
    xlim([0 60])
    ylim([0 ylimsup(m)]);
    ylabel('Weekly new cases','FontSize',10)
    set(gca,'FontSize',8,'XColor','k','YColor','k','TickLength',[0.025 0.025])
    print(h,strcat('auxinset_no_',num2str(m)),'-dpdf')
   
end


    setup.R0max = [R0v1 R0max_swipe(I)]';
    setup.Tin   = [-5 500];
    [Par,Fun]   = generador_parametros(setup);
    odefun = @(t,y) dydt(t,y,Par,Fun);
    generador_cond_inicial(Par)
    [t,y] = ...
        ode45(@(t,y) odefun(t,y), tspan, generador_cond_inicial(Par), options);
    [Nobs_com, Nobs_POE, N_com, N_POE] = Nobs(t,y,Par,Fun);
    %% Accumulating and discretising weekly cases
    
    Nobs_com_weekly = nan(k,floor((setup.tf-setup.ti)/7) + 1);
    Nobs_POE_weekly = nan(k,floor((setup.tf-setup.ti)/7) + 1);
    T = floor(setup.ti/7):floor(setup.tf/7);
    for i = 2:floor((setup.tf-setup.ti)/7) + 1
        for j = 1:k
            idx = t/7 >= T(i-1) & t/7 <= T(i);
            Nobs_com_weekly(j,i) = floor(trapz(t(idx),Nobs_com(j,idx)));
            Nobs_POE_weekly(j,i) = floor(trapz(t(idx),Nobs_POE(j,idx)));
        end
    end
    
    
    h = figure('units','centimeters','position',[3,3,9,3]);
    plot(t/7,Nobs_com(1,:),'LineWidth',2)
    hold on
    xlim([10 50])
    ylim([0 300]);
    ylabel('Weekly new cases','FontSize',10)
    set(gca,'FontSize',8,'XColor','k','YColor','k','TickLength',[0.025 0.025])
    print(h,strcat('auxinset_no_',num2str(4)),'-dpdf')

