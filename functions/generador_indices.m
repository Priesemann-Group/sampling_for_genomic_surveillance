function ID = generador_indices(k)

ID_inner = logical(repmat(eye(5),k,1));
ID = false(5*k+2,5+2); ID(1,1) = true; ID(end,end) = true;
ID(2:end-1,2:end-1) = ID_inner;