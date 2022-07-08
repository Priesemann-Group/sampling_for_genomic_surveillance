function fi_hat = GS(Nobsx_week,kx,m)

NTx_week = sum(Nobsx_week);
if kx < NTx_week
    rng(m)
    indGS = randperm(NTx_week,kx);
    fi_hat = zeros(size(Nobsx_week));
    idx = ~(Nobsx_week==0);
    Nfil = Nobsx_week(idx);
    Vk_aux = zeros(1,sum(idx));
    for j = 1:sum(idx)
        Vk_aux(j) = sum(indGS<=sum(Nfil(1:j)));
    end
    Vk_aux = Vk_aux - [0 Vk_aux(1:end-1)];
    fi_hat(idx) = Vk_aux;
    fi_hat = fi_hat/sum(fi_hat);
elseif NTx_week > 0
    fi_hat = Nobsx_week/sum(Nobsx_week);
else
    fi_hat = zeros(size(Nobsx_week));
end
    
    