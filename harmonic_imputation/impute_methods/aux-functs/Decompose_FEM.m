function [IMF,theta,phi,dtheta]=Decompose_FEM(f,theta_ini)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%{

% This is a function do seperate one IMF, env*cos(theta), from given data f.


% INPUT:
% f: original data
% theta_ini: the initial value of the phase function theta. It must have
% same size as f and evaluated over uniform grids in [0,1] since in this program
% we assume the span of data in time is [0,1].
%
% OUTPUT:
% IMF: the intrinsic mode function extracted by our method;
% theta: phase function;
% phi: phase constant since phase function theta is enforced to be 0 at the
% original point t=0;
% dtheta: instantaneous frequency which is the derivative of theta.




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% initialization of some working variables

x=linspace(0,1,length(f))';
N=length(x);
c=floor(log(N)/log(2));
N_t=2^c;

% initial value of theta

theta=theta_ini;
dtheta=diff_center(theta);

M=40;%max(round((theta(end)-theta(1))/(2*pi)),30);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

% alpha is the parameter to control the smoothness of theta
% first alpha is small to limit theta to a smooth space
% then this limitation is relaxed gradually.

for alpha=1:c
    if (2^(alpha-1)>(theta(end)-theta(1))/(4*pi))
        break;
    end
    mesh_b=linspace(0,1,2^(alpha-1)+1);
    M_b=basis_FEM_linear(mesh_b,x);
    for step=1:10 
        
        % interpolate the data to the theta space
        h_t=(theta(end)-theta(1))/N_t;
        theta_fit=[theta(1):h_t:theta(end)]';
        f_fit=spline(theta,f,theta_fit);
        
        
        
        % extract the IMF and a, b by FFT in theta space
        [a1_c,a1_s,IMF_fit]=decompose_env_FEM(f_fit,theta_fit,theta_fit/theta(end),N_t);
        


        theta_old=theta;
        dtheta_old=dtheta;

        % interpolate the results to time space
        a1_ct=spline(theta_fit,a1_c,theta);
        a1_st=spline(theta_fit,a1_s,theta);
        IMF=spline(theta_fit,IMF_fit,theta);

        % calculate the envelop
        
        env2=(a1_ct.^2+a1_st.^2);
        env=sqrt(env2);

        % calculate the derivative by centeral difference
        
        da1_ct=diff_center(a1_ct);
        da1_st=diff_center(a1_st);

        % calculate the change of theta
        
        env2=max(env2,max(env2)/100);
        dtheta_add=(a1_ct.*da1_st-a1_st.*da1_ct)./env2;


        % apply low-pass filter on the update of theta
        % the smoothness is controlled by the parameter alpha
        
        %dtheta_add_s = dtheta_add;
        
        dtheta_add_s=smooth_IF_FEM(dtheta_add,M_b,mesh_b,x,alpha);
        

        % calculate the parameter that make the frequency positive
        id=dtheta_add_s>0;
        if(sum(id)==0)
            lambda=1;
        else
            lambda=min(min(dtheta(id)./dtheta_add_s(id))/2,1);
        end

        % update the frequency
        dtheta=2*pi*round((theta(end)-theta(1))/(2*pi))/(theta(end)-theta(1))*dtheta-lambda*dtheta_add_s(1:N_t);

        % reconstruct the phase function from the frequency
        theta=integral_fft_sym(dtheta);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        % print the error

        error=norm(dtheta-dtheta_old,2)/norm(dtheta,2);
        %fprintf('al=%1.3f, step=%d, error=%1.2e, res=%e, lambda=%1.2f\n',alpha, step, error, norm(f-IMF,2),lambda);


        if (error<1e-2)
            break;
        end
    end
end

%determine the phase constant
            
v_theta=theta-theta_old;

rc=a1_ct.*cos(v_theta)-a1_st.*sin(v_theta);
rs=-a1_st.*cos(v_theta)-a1_ct.*sin(v_theta);

cp=sum(env.*rc)/(sum(env.^2));
sp=sum(env.*rs)/(sum(env.^2));

phi=atan(sp/cp);
if (cp<0)
    phi=phi+pi;
end