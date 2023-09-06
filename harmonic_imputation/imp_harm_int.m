function [s_imp,trend_imp,s_int,trend_int,bestM, A_int, phi_int,TFData] = imp_harm_int(s,varargin)
%%
% Compute the harmonic decomposition of signal x.
% Inputs: 
%         s: signal under analysis
%         Additional arguments:
%           Int: String variable. Set to "On" to perform harmonic
%           interpolation
%           Best: String variable. Set to "On" to choose the best
%           performing imputation method to carry over to the harmonic
%           decomposition step. Needs to specify ground truth.
%           True: Ground truth used to measure the error when comparing
%           imputation methods.
%           ImpM: String of Cell Array. Indicates imputation method to
%           compute.
%           params_mss: Struct. Indicates the parameters for the missing
%           interval detection algorithm. See missing_ints.m for details.
%           params_imp: Struct. Indicates the parameters for the initial 
%           imputation step algorithm. See impute.m for details.
%           params_decomp: Struct. Indicates the parameters for the 
%           harmonic decomposition step. See harm_decomp.m for details.
%           IntM: String or Cell Array. Choose which interpolation methods
%           to use at the harmonic interpolation step.
%
% Outputs:
%         s_imp: Initial Imputation result for each method.
%         s_int: Harmonic interpolation result for each method and
%         interpolator.
%

p = inputParser;
addOptional(p,'Int','Off');
addOptional(p,'Best','Off');
addOptional(p,'True',0);
addOptional(p,'ImpM','tlm');
addOptional(p,'params_mss',struct());
addOptional(p,'params_imp',struct());
addOptional(p,'params_decomp',struct());
addOptional(p,'IntM',{'pchip'});
addOptional(p,'params_missints',struct());
parse(p,varargin{:})

if isfield(p.Results.params_mss,'st') && isfield(p.Results.params_mss,'L')
    st = p.Results.params_mss.st;
    L =  p.Results.params_mss.L;
else
    if isfield(p.Results.params_mss,'d')
        d =  p.Results.params_mss.d;
    else
        d = 0.01*length(s);
    end
    if isfield(p.Results.params_mss,'t')
        t =  p.Results.params_mss.t;
    else
        t = 0;
    end
    [st,L] = missing_ints(s,p.Results.params_missints);
end

s_imp = impute(s,st,L,p.Results.ImpM,p.Results.params_imp);

if strcmpi(p.Results.Best,'on')
    mse_imp = zeros(1,size(s_imp));

    for k=1:size(s_imp,1)
        si = s_imp(k,:)';
        [mse_imp(k)] = compute_errors(p.Results.t,si,st,L);
    end
    [~,bestM] = min(mse_imp);

    s_imp = s_imp(bestM,:);

else
    bestM = 0;
end

if  strcmpi(p.Results.Int,'on')
    K = size(s_imp,1);
    M = size(p.Results.IntM,2);

    s_int = zeros(M*K,length(s));
    trend_int = zeros(M*K,length(s));
    A_int = cell(1,M*K);
    phi_int = cell(1,M*K);
    for k=1:size(s_imp,1)
        si = s_imp(k,:)';
        [AD,phiD,trend_imp,TFData] = harm_decomp(si, p.Results.params_decomp);
        for m=1:M
            intmm =  p.Results.IntM{m};

            [s_int((k-1)*M+m,:),trend_int((k-1)*M+m,:),...
             A_int{(k-1)*M+m},phi_int{(k-1)*M+m}] = harm_int(AD,phiD,st,L,intmm,si,trend_imp);
        end
    end
end

end