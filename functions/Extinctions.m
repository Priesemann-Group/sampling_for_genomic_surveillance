function [value,isterminal,direction] = Extinctions(t,y,Par)
% Integrador de las ecuaciones diferenciales del sistema TI multivariantes
% de SARS-CoV-2;

%% Definición de variables de estado

EQ_i    = y(Par.ID(:,2),:);
EH_i    = y(Par.ID(:,3),:);
IHs_i   = y(Par.ID(:,4),:);
IHa_i   = y(Par.ID(:,5),:);
IQ_i    = y(Par.ID(:,6),:);

Cases_i = EQ_i + EH_i + IHs_i + IHa_i + IQ_i;
extinctions = Cases_i-Par.cutOff*ones(size(Cases_i));

% The value that we want to be zero
% Seems redundant, but I keep this line as a memo: New events can be added
value = extinctions;%[extinctions];

% Halt integration if some value reach zero
isterminal = ones(size(value));
% The zero can be approached:
%     0: from both directions
%     1: when function is increasing
%    -1: when function is decreasing
direction = -1*ones(size(value));

end
