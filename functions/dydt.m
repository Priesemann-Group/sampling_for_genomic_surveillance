function F = dydt(t,y,Par,Fun)
% Integrador de las ecuaciones diferenciales del sistema TI multivariantes
% de SARS-CoV-2;

%% Definición de variables de estado

S       = y(Par.ID(:,1),:);
EQ_i    = y(Par.ID(:,2),:);
EH_i    = y(Par.ID(:,3),:);
IHs_i   = y(Par.ID(:,4),:);
IHa_i   = y(Par.ID(:,5),:);
IQ_i    = y(Par.ID(:,6),:);

%% Obtención de varibles temporales

Rt_i    = Fun.Rt(t,Par);
Phi_i   = Fun.Phi(t,Par);

%% Derivadas

dSdt    = -Par.gamma*(S/Par.M)*dot(Rt_i,(IHs_i+IHa_i)) ...
    -Par.gamma*(S/Par.M)*(Par.epsilon+Par.nu)*dot(Rt_i,IQ_i) - (S/Par.M)*sum(Phi_i);
dEQ_i   = Par.gamma*(S/Par.M)*Par.nu*(Rt_i.*IQ_i) - Par.rho * EQ_i;
dEH_i   = Par.gamma*(S/Par.M)*Par.epsilon*(Rt_i.*IQ_i) ...
    + Par.gamma*(S/Par.M)*(Rt_i.*(IHs_i+IHa_i)) - Par.rho * EH_i;
dIHs_i  = Par.rho*Par.xi.*EH_i - Par.gamma*IHs_i ...
    -(Par.lambda_s + Par.lambda_r)*Par.eta .* IHs_i ...
    + (Par.xi.*(1-Par.eta)).* Phi_i *(S/Par.M);
dIHa_i  = Par.rho*(1-Par.xi).*EH_i - Par.gamma*IHa_i ...
    -Par.lambda_r*Par.eta .* IHa_i ...
    + ((1-Par.xi).*(1-Par.eta)).* Phi_i *(S/Par.M);
dIQ_i   = Par.rho*EQ_i - Par.gamma*IQ_i ...
    + (Par.lambda_s + Par.lambda_r)*Par.eta .* IHs_i ...
    + Par.lambda_r*Par.eta .* IHa_i + Par.eta .* Phi_i * S/Par.M;
dRdt    = Par.gamma*sum(IQ_i + IHs_i + IHa_i);

%% Reordering;

F = zeros(length(Par.ID(:,1)),1);
F(Par.ID(:,1)) = dSdt;
F(Par.ID(:,2)) = dEQ_i;
F(Par.ID(:,3)) = dEH_i;
F(Par.ID(:,4)) = dIHs_i;
F(Par.ID(:,5)) = dIHa_i;
F(Par.ID(:,6)) = dIQ_i;
F(Par.ID(:,7)) = dRdt;
