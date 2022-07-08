function N_seq_obs = cuenta_secuencias(Nobs,ID,k)
% función que cuenta cuántas muestras de las k_seq secuenciadas en un
% tiempo determinado corresponden a cada una de las k variantes. Nobs es un
% vector con k entradas, representando los casos observados para cada
% variante luego de un delay de testeo. 

Nobs = floor(Nobs);
M = Nobs'*triu(ones(k))';
ni = zeros(k,1);

for i = 1:k
    ni(i) = sum(ID<=M(i));
end

N_seq_obs = ni'.*inv(triu(ones(k)));