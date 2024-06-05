function ALP = sample_haf(t,subs,n)

if nargin<3
    n = 1;
end

if nargin<2
    subs = {'Poly','Trigo','Tanh','Bump'};
end

N = length(t);

ALP = zeros(n,N);
for k=1:n
    i = randi(size(subs));
    s = subs{i};

    switch s
        case 'Const'
            alp = 0.8*rand(1,1)*ones(1,length(t));
        case 'Poly'
            ord = randi(4);
            coefs = rand(1,4);
            coefs = coefs/sum(coefs)-0.1;
            alp = 0;
            for k=1:ord
                alp = alp + coefs(k)*(t/t(end)).^k;
            end

        case 'Trigo'
            m = 0.1 + 0.7*rand(1,1);
            a = 0.1 + (1-m-0.1)*rand(1,1);
            f = randi(6);
            if randi(2)<1
                alp = m + a*cos(2*pi*f*t);
            else
                alp = m + a*sin(2*pi*f*t);
            end

        case 'Tanh'
            t0 = 0.25 + 0.5*rand(1);
            p  = 0.05 + 0.1*rand(1);
            m = rand(1);
            alp = 0.5 + p + 0.25*tanh(sqrt(N)*m*(t-t0));

        case 'Bump'
            sig = rand(1);
            t0 = 0.25 + 0.5*rand(1);
            m = 0.1 + 0.6*rand(1,1);
            alp = m + 0.25*exp(-((t-t0)./sig).^2);
        otherwise
            alp = zeros(1,N);
    end

    ALP(k,:) = alp;
end

lambd = rand(1,n);
lambd = lambd/sum(lambd);

end