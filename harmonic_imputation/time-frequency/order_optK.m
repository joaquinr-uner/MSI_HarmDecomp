function [r_opt,Crits] = order_optK(s,r_max,A,phi,criteria,params,F)
%% 
% Optimum number of harmonics estimation by trigonometric model regression
% selection criteria for multicomponent signals 
% Inputs:
%         s: signal under analysis
%         r_max: maximum admissible number of harmonics.
%         A: instantaneous amplitude (IA) of s.
%         phi: instantaneous phase (IP) of s.
%         criteria: model selection criteria to be used. Select from
%                  'GCV': Generalized Cross-Validation.
%                  'Rl': Unbiased Risk.
%                  'Wang': Wang Criterion.
%                  'Kavalieris': Kavalieris-Hannan Criterion.
%         params: model selection criteria parameters.
%                  'Rl': Noise estimator sigma.
%                  'Wang': penalization parameter c.
%                  'Kavalieris': penalization parameter H.
RSS = zeros(1,prod(r_max));

[K,N] = size(A);
if size(s,1) == 1
    s = s';
end
if nargin<7
    F = [];
end

if nargin<6
   params = struct('None',[]); 
    
end
if nargin<5
    criteria = {'GCV'};
else
    if ~iscell(criteria)
        criteria = {criteria};
    end
end

for i=1:length(criteria)
    cri=criteria{i};
    switch cri
        case 'GCV'
            GCV = zeros(1,prod(r_max));
        case 'Rl'
            if length(criteria) == 1 && ~isa(params,'struct')
                sigma = params;
            else
                if isfield(params,'sigma')
                    sigma = params.sigma;
                else
                    if ~isempty(F)
                        sigma = sqrt(2)*median(abs(real(F(:))))/0.6745;
                    else
                        fprintf('STFT Needed')
                    end
                end
            end
            Rl = zeros(1,prod(r_max));
        case 'Wang'
            if length(criteria)==1 && ~isa(params,'struct')
                rc = params;
            else
                if isfield(params,'c')
                    rc = params.c;
                else
                    rc = 2.1;
                end
            end
            Wn = zeros(prod(r_max),length(rc));
        case 'Kavalieris'
            if length(criteria)==1 && ~isa(params,'struct')
                H = params;
            else
                if isfield(params,'H')
                    H = params.H;
                else
                    H = round(log(N)^2);
                end
            end
            Kv = zeros(prod(r_max),H);
        otherwise
            fprintf('%s is not a valid selection criteria \n',cri)
    end
end
if ~isa(criteria,'cell')
    criteria = {criteria};
end

C = construct_dct(A,phi,r_max);

rcell = cell(1,K);
for k=1:K
    rcell{k} = 1:r_max(k); 

end

vecr = combvec(rcell{1,:});

for r=1:size(vecr,2)
    Cr = construct_dct(A,phi,vecr(:,r));
    v = ((Cr'*Cr)\Cr')*s;

    s_est = Cr*v;
    E = (s-s_est);
    RSS(r) = sum(E.^2);
    for i=1:length(criteria)
        cri = criteria{i};
        switch cri
            case 'GCV'
                GCV(r) = N*RSS(r)/(N - 2*sum(vecr(:,r)) - 1)^2;
            case 'Rl'
                Rl(r) = 1/N*RSS(r) + 2*(sigma)^2*(2*sum(vecr(:,r)) + 1)/N;
            case 'Wang'
                E = (s-s_est);
                RSS(r) = sum(E.^2);
                for ci=1:length(rc)
                    Wn(r,ci) = log10(1/N*RSS(r)) + rc(ci)*sum(vecr(:,r))*log10(N)/N;
                end
            case 'Kavalieris'
            for h=1:H
                [~,var] = lpc(E,h);
                Kv(r,h) = log(var) + (5*sum(vecr(:,r)) + h)*log10(N)/N;
            end
        end  
    end
end

r_opt = zeros(K,length(criteria));
Crits = struct();
for i=1:length(criteria)
    cri = criteria{i};
    switch cri
        case 'GCV'
            [~,aux] = min(GCV);
            r_opt(:,i) = vecr(:,aux);
            Crits.GCV = GCV;
        case 'Rl'
            [~,aux] = min(Rl);
            r_opt(:,i) = vecr(:,aux);
            Crits.Rl = Rl;
        case 'Wang'
            [~,aux] = min(Wn(:));
            [bc,~] =  ind2sub(size(Wn),aux);
            r_opt(:,i) = vecr(:,bc);
            Crits.Wn = Wn;
        case 'Kavalieris'
            [~,aux] = min(Kv(:));
            [bc,~] =  ind2sub(size(Kv),aux);
            r_opt(:,i) = vecr(:,bc);
            Crits.Kv = Kv;
    end
end
