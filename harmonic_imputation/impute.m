function [x_imp,t_imp] = impute(x,st,L,method,params,verbose)
%%
% Wrapper function for imputation methods. Can be used to run
% multiple imputation methods at once. Returns a matrix with each row
% containing an imputed signal.
%
% Inputs:
%         x: signal with missing data.
%         st: initial indexes of missing data intervals.
%         L: lengths of missing data intervals.
%         method: either a cell of strings with the names of the
%                 built-in methods or a function handle to a custom
%                 imputation method. The custom imputation method has
%                 to follow the standard definition:
%                   function [s_imp] = method(s,st,L,params)
%         params: optional parameter struct or cell of structs.
% Outputs:
%       x_imp: signals with imputed values.
%       t_imp: computation time of the imputation methods.

if nargin<6
    verbose = false;
end

if nargin<5
    params = '0';
end
if nargin<4
    tic
    x_imp = impute_tlm(x,st,L);
    t_imp = toc;
else
    if iscell(method)
        if isa(params,'struct')
            params = {params};
        end

        K = size(method,2);
        x_imp = zeros(K,length(x));
        t_imp = zeros(K,1);
        for k=1:K
            switch lower(method{k})
                case 'tlm'
                    if nargin<5
                        tic
                        x_imp(k,:) = impute_tlm(x,st,L);
                        t_imp(k) = toc;
                        if verbose
                            fprintf(['Imputation completed by TLM method. Time: ' num2str(t_imp(k)) '\n'])
                        end
                    else
                        tic
                        x_imp(k,:) = impute_tlm(x,st,L,params{k});
                        t_imp(k) = toc;
                        if verbose
                            fprintf(['Imputation completed by TLM method. Time: ' num2str(t_imp(k)) '\n'])
                        end
                    end
                case 'lse'
                    if nargin<5
                        tic
                        x_imp(k,:) = impute_lse(x,st,L);
                        t_imp(k) = toc;
                        if verbose
                            fprintf(['Imputation completed by LSE method. Time: ' num2str(t_imp(k)) '\n'])
                        end
                    else
                        tic
                        x_imp(k,:) = impute_lse(x,st,L,params{k});
                        t_imp(k) = toc;
                        if verbose
                            fprintf(['Imputation completed by LSE method. Time: ' num2str(t_imp(k)) '\n'])
                        end
                    end
                case 'dmd'
                    if nargin<5
                        tic
                        x_imp(k,:) = impute_dmd(x,st,L);
                        t_imp(k) = toc;
                        if verbose
                            fprintf(['Imputation completed by DMD method. Time: ' num2str(t_imp(k)) '\n'])
                        end
                    else
                        tic
                        x_imp(k,:) = impute_dmd(x,st,L,params{k});
                        t_imp(k) = toc;
                        if verbose
                            fprintf(['Imputation completed by DMD method. Time: ' num2str(t_imp(k)) '\n'])
                        end
                    end
                case 'gpr'
                    if nargin<5
                        tic
                        x_imp(k,:) = impute_gpr(x,st,L);
                        t_imp(k) = toc;
                        if verbose
                            fprintf(['Imputation completed by GPR method. Time: ' num2str(t_imp(k)) '\n'])
                        end
                    else
                        tic
                        x_imp(k,:) = impute_gpr(x,st,L,params{k});
                        t_imp(k) = toc;
                        if verbose
                            fprintf(['Imputation completed by GPR method. Time: ' num2str(t_imp(k)) '\n'])
                        end
                    end
                case 'arimaf'
                    if nargin<5
                        tic
                        x_imp(k,:) = impute_arimaf(x,st,L);
                        t_imp(k) = toc;
                        if verbose
                            fprintf(['Imputation completed by ARIMAF method. Time: ' num2str(t_imp(k)) '\n'])
                        end
                    else
                        tic
                        x_imp(k,:) = impute_arimaf(x,st,L,params{k});
                        t_imp(k) = toc;
                        if verbose
                            fprintf(['Imputation completed by ARIMAF method. Time: ' num2str(t_imp(k)) '\n'])
                        end
                    end
                case 'arimab'
                    if nargin<5
                        tic
                        x_imp(k,:) = impute_arimab(x,st,L);
                        t_imp(k) = toc;
                        if verbose
                            fprintf(['Imputation completed by ARIMAB method. Time: ' num2str(t_imp(k)) '\n'])
                        end
                    else
                        tic
                        x_imp(k,:) = impute_arimab(x,st,L,params{k});
                        t_imp(k) = toc;
                        if verbose
                            fprintf(['Imputation completed by ARIMAB method. Time: ' num2str(t_imp(k)) '\n'])
                        end
                    end
                case 'ddtfa'
                    if nargin<5
                        tic
                        x_imp(k,:) = impute_ddtfa(x,st,L);
                        t_imp(k) = toc;
                        if verbose
                            fprintf(['Imputation completed by DDTFA method. Time: ' num2str(t_imp(k)) '\n'])
                        end
                    else
                        tic
                        x_imp(k,:) = impute_ddtfa(x,st,L,params{k});
                        t_imp(k) = toc;
                        if verbose
                            fprintf(['Imputation completed by DDTFA method. Time: ' num2str(t_imp(k)) '\n'])
                        end
                    end
                case 'tbats'
                    if nargin<5
                        tic
                        x_imp(k,:) = impute_tbats(x,st,L);
                        t_imp(k) = toc;
                        if verbose
                            fprintf(['Imputation completed by TBATS method. Time: ' num2str(t_imp(k)) '\n'])
                        end
                    else
                        tic
                        x_imp(k,:) = impute_tbats(x,st,L,params{k});
                        t_imp(k) = toc;
                        if verbose
                            fprintf(['Imputation completed by TBATS method. Time: ' num2str(t_imp(k)) '\n'])
                        end
                    end
                case 'edmd'
                    if nargin<5
                        tic
                        x_imp(k,:) = impute_edmd(x,st,L);
                        t_imp(k) = toc;
                        if verbose
                            fprintf(['Imputation completed by EDMD method. Time: ' num2str(t_imp(k)) '\n'])
                        end
                    else
                        tic
                        x_imp(k,:) = impute_edmd(x,st,L,params{k});
                        t_imp(k) = toc;
                        if verbose
                            fprintf(['Imputation completed by EDMD method. Time: ' num2str(t_imp(k)) '\n'])
                        end
                    end
                case 'lsw'
                    if nargin<5
                        tic
                        x_imp(k,:) = impute_lsw(x,st,L);
                        t_imp(k) = toc;
                        if verbose
                            fprintf(['Imputation completed by LSW method. Time: ' num2str(t_imp(k)) '\n'])
                        end
                    else
                        tic
                        x_imp(k,:) = impute_lsw(x,st,L,params{k});
                        t_imp(k) = toc;
                        if verbose
                            fprintf(['Imputation completed by LSW method. Time: ' num2str(t_imp(k)) '\n'])
                        end
                    end
                otherwise
                    error(['Unknown imputation method: ' method{k}])
            end

        end
    else
        if nargin<5
            x_imp = method(x,st,L);
        else
            x_imp = method(x,st,L,params);
        end
    end
end
