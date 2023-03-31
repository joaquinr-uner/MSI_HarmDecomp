function [Tmse,mse,Tmae,mae,Trmse,rmse] = compute_errors(t,x,st,L)
    K = length(L);

    mse = zeros(1,K);
    mae = zeros(1,K);
    rmse = zeros(1,K);
    int = [];
    for k=1:K
        Lk = L(k);
        stk = st(k);
        edk = stk + Lk - 1;
        intk = stk:edk;
        mse(k) = norm(t(intk)-x(intk))^2/Lk;
        mae(k) = norm(t(intk)-x(intk),1)/Lk;
        rmse(k) = sqrt(mse(k));
        int = [int intk];


    end
    Tmse = norm(t(int)-x(int))^2/sum(L);
    Tmae = norm(t(int)-x(int),1)/sum(L);
    Trmse = sqrt(Tmse);




end

