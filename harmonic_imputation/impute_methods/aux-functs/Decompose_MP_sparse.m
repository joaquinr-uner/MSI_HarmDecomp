function [a0,IMF,theta_fine,env,dtheta]=Decompose_MP_sparse(f,theta_ini,xs,xf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%{

% This is a function to extract IMF, env*cos(theta), from incomplete or sparse data f.
% theta_ini is the initial value of the phase function theta.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% initialization of some working variables

p2=2*pi;


theta_fine=theta_ini;
theta=spline(xf,theta_fine,xs);
dtheta=diff_center(theta_fine);

M=10;%max(round((theta(end)-theta(1))/(2*pi)),30);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

for alpha=1/floor(M):1/floor(M):1/5
    theta_ori=theta;
    converge=0;
    for beta=1:1
        theta=theta_ori;
        for step=1:10
                      
            [a1_ct,a1_st,a0,IMF]=decompose_sparse(f,theta,theta_fine,1/2*2^(1-beta));
            
            dtheta_old=dtheta;

            
            env2=(a1_ct.^2+a1_st.^2);
            env=sqrt(env2);

            da1_ct=diff_center(a1_ct);
            da1_st=diff_center(a1_st);


            env2=max(env2,0.01);
            dtheta_add=(a1_ct.*da1_st-a1_st.*da1_ct)./env2;



            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
            % some processing to make sure d theta/dt>0
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            
            dtheta_add_s=smooth_theta(dtheta_add(1:end-1),theta_fine,alpha);
            dtheta_add_s=[dtheta_add_s;dtheta_add_s(1)];

            % calculate the parameter that make the frequency positive
            id=dtheta_add_s>0;
            if(sum(id)==0)
                lambda=1;
            else
                lambda=min(min(dtheta(id)./dtheta_add_s(id))/2,1);
            end
            
            % update the frequency
            dtheta=dtheta-lambda*dtheta_add_s;



            % reconstruct the phase function from the frequency
            theta_fine=integral_fft(dtheta(1:end-1));%
            theta=spline(xf,theta_fine,xs);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
            % plot the frequency we get in each time step
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
            error=norm(dtheta-dtheta_old,2)/norm(dtheta,2);



            if (min(dtheta)<0 || error<1e-2)
                converge=1;
                break;
            end
        end
        if converge>0
            break;
        end
    end
end