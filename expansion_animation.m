function ret = expansion_animation(t,Y,ds,rho,c_t,c_a,R_x_inds,R_y_inds,dR_x_inds,dR_y_inds,P_ext,R0,fig,fps,tspan,filename)
    n = size(t,1)
    h = figure(fig)
    frames = 1+tspan(1)*fps:tspan(2)*fps
    axis tight manual % this ensures that getframe() returns a consistent size
    for i = frames
        clf(fig)
        hold on
        for j = 1:size(Y,3)
            y = Y(:,:,j);
            [~,t_ind] = min(abs(t-i/fps))
            [th,k,dsds,k_dot,dsds_dot] = get_dshapedt([y(t_ind,R_x_inds).',y(t_ind,R_y_inds).'],[y(t_ind,dR_x_inds).',y(t_ind,dR_y_inds).'],ds);
            plot(y(t_ind,R_x_inds),y(t_ind,R_y_inds))
            F_bend = [k_t(0)/(rho)*k(:,1).';k_t(0)/(rho)*k(:,2).']*ds^2;
            F_ax = [k_a(0)/(rho)*(dsds(:,1)).';k_a(0)/(rho)*(dsds(:,2)).']*ds^2;
            C_bend = [c_t*k_dot(:,1).';c_t*k_dot(:,2).']*ds^2;
            C_ax = [c_a*dsds_dot(:,1).';c_a*dsds_dot(:,2).']*ds^2;
            quiver(y(t_ind,R_x_inds),y(t_ind,R_y_inds),F_bend(1,:)- F_ax(1,:),F_bend(2,:)- F_ax(2,:),'red','AutoScale','off')
            quiver(y(t_ind,R_x_inds),y(t_ind,R_y_inds), C_bend(1,:)- C_ax(1,:),C_bend(2,:)- C_ax(2,:),'magenta','AutoScale','off')
            quiver(y(t_ind,R_x_inds),y(t_ind,R_y_inds),-P_ext(t(t_ind)).'.*abs(th(:,2)).'*ds^2.*y(t_ind,R_x_inds),-P_ext(t(t_ind)).*abs(th(:,1)).'*ds^2.*y(t_ind,R_y_inds),'cyan','AutoScale','off')
            title(strcat("Expansion Simulation, no External Forces, t=",num2str(t(t_ind))))
            ylim([-2*R0 2*R0])
            xlim([-2*R0 2*R0])
            drawnow
            t(t_ind)
        end
        hold off
    end
    %{
    for i = 1:n-1
        i
            [A,map] = rgb2ind(im{i},256);
            if i == 1
                imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0);
            else
                imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0);
            end
     end
    %}
end