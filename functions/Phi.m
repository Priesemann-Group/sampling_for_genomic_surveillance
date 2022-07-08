function Phi_i = Phi(t,Par)
Phi_i = zeros(Par.k,length(t));
idx = true(Par.k,length(t));
min_time = Par.min_time;
Phibase = Par.Phibase;
Tin = [Par.Tin Par.tspan(end) Par.tspan(end)];
Tin(1) = Tin(1) - min_time;

for i = 1:Par.k
    Phi_i(i,:) = Par.Phimax(i)*gampdf(t-Par.Tin(i),Par.a(i),Par.b(i));
    idx(i,:) = t >= Tin(i) + min_time & t <= Tin(i+1) + 2*min_time;
end

%% fix para incluír siempre la última variante en un influx razonable

Phi_i = max(Phi_i,Phibase * idx./sum(idx));

