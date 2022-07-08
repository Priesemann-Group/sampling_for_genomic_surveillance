function [Nobs_com, Nobs_POE, N_com, N_POE] = Nobs(t,y,Par,Fun)
% Función que entrega los casos nuevos observados en aeropuerto y en la
% comunidad. Esta cantidad se calcula apartando la parte positiva de la
% derivada en y(t) y convolucionándola con un kernel de delay (empírico o
% gamma)

%% Definición de variables

y = y';
S       = y(Par.ID(:,1),:);
IHs_i   = y(Par.ID(:,4),:);
IHa_i   = y(Par.ID(:,5),:);

%% Obtención de varibles temporales
Phi_i   = Fun.Phi(t,Par);

%% Sampling at POE

N_POE = Par.eta .* (Phi_i .* S/Par.M);
Nobs_POE = N_POE;
for i = 1:Par.k
    Nobs_POE(i,:) = EstimDelay(Nobs_POE(i,:),4,1,0.95);
end
%% Derivadas

N_com = (Par.lambda_s + Par.lambda_r)*Par.eta .* IHs_i + Par.lambda_r*Par.eta .* IHa_i;
Nobs_com = N_com;
for i = 1:Par.k
    Nobs_com(i,:) = EstimDelay(Nobs_com(i,:),4,1,0.95);
end