function x0 = generador_cond_inicial(Par)

x0 = zeros(length(Par.ID(:,1)),1);

x0(Par.ID(:,2)) = Par.isVarThere.*Par.I0*Par.nu;              % EQ
x0(Par.ID(:,3)) = Par.isVarThere.*Par.I0*(1+Par.epsilon);                                   % EH
x0(Par.ID(:,4)) = Par.isVarThere.*Par.I0*(Par.nu+Par.epsilon);              % IQ
x0(Par.ID(:,5)) = Par.isVarThere.*Par.I0.*Par.xi;                           % IHs
x0(Par.ID(:,6)) = Par.isVarThere.*Par.I0.*(1-Par.xi);                       % IHa
% x0(Par.ID(:,7)) = Par.isVarThere.*Par.I0*(Par.nu+Par.epsilon);
% % R = 0
x0(Par.ID(:,1)) = Par.M-sum(x0(2:end));
