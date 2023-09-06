function [st,ed,Ll] = rand_missing_ints(L,N,Ni)

    Lmin = round(L/4);
    Lmax = round(L/2);
    Ll = floor(randfixedsum(Ni,1,L,Lmin,Lmax)');
    st = zeros(Ni,1);
    ed = zeros(Ni,1);
    for i=1:Ni
        st(i) = round((i/(Ni+1)-0.05)*N) + randi(0.1*N);
        ed(i) = st(i) + Ll(i);
    end