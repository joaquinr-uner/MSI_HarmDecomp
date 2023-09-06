function [st,LEN,Ni] = missing_ints(s, params)
%%
% Find missing data intervals in signal s, according to the missing data
% classification proposed in Ruiz, J, et al. "Enhanching Missing Data..."
% Inputs:
%         s: signal with missing data.
%         params: struct variable with missing data algorithm parameters.
%         Valid fields listed below:
%            c: missing data type to be considered. Valid values are:
%                   'x': x-missing case (disconnection).
%                   'y': y-missing case (oversaturation).
%                   'c': corruptive case.
%
%            t: threshold (for x-missing case).
%            d: minimun missing data interval (for x-missing and y-missing cases).
%            sqi: signal quality index (for corruptive case).
%            q0: missing admmisible quality value (for corruptive case).
% Outputs:
%         st: initial indexes of missing data intervals.
%         LEN: lengths of missing data intervals.
%         Ni: number of distinct missing data intervals.
if nargin<2 && isempty(~isfield(params,'c'))
    params = struct('c','x','d',5,'t',0);
end
if isfield(params,'c')
    c = params.c;
else
    c = 'x';
end
if isfield(params,'d')
    d = params.d;
else
    %b = round(sqrt(log(80)*sigma/pi^2)*length(x)) + 1;
    d = 0.01*length(s);
end
if isfield(params,'t')
    t = params.t;
else
    t = 0;
end
if isfield(params,'sqi')
    sqi = params.sqi;
else
    sqi = Inf*ones(size(s));
end
if isfield(params,'q0')
    q0 = params.q0;
else
    q0 = Inf;
end

switch c
    case 'x'
        tmp = abs(s)<=t;
        [M, st] = regexp(sprintf('%i', tmp'), '1+', 'match') ;

        p = 0;
        for i = 1:size(M,2)
            if length(M{i})>d
                p = p + 1;
                LEN(p) = length(M{i});
                aux(p) = st(i);
            end
        end
        st = aux;
        Ni = length(LEN);
    case 'y'
        tmp = s==max(s);
        [M, st] = regexp(sprintf('%i', tmp'), '1+', 'match') ;

        p = 0;
        for i = 1:size(M,2)
            if length(M{i})>d
                p = p + 1;
                LEN(p) = length(M{i});
                aux(p) = st(i);
            end
        end
        st = aux;
        Ni = length(LEN);
    case 'c'
        tmp = s(sqi<q0);
        [M, st] = regexp(sprintf('%i', tmp'), '1+', 'match') ;

        p = 0;
        for i = 1:size(M,2)
            if length(M{i})>d
                p = p + 1;
                LEN(p) = length(M{i});
                aux(p) = st(i);
            end
        end
        st = aux;
        Ni = length(LEN);
    otherwise

        error('Invalid missing data type')
end
end