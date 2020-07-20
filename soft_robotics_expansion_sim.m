S = 10 %mm
ds = .01 %mm
n = round(S/ds)
rho = .70

inds = [1:n*8].';
R_x_inds = inds(mod(1:size(inds,1),8)==1)
th_x_inds = inds(mod(1:size(inds,1),8)==2);
k_x_inds = inds(mod(1:size(inds,1),8)==3);
dR_x_inds = inds(mod(1:size(inds,1),8)==4);
R_y_inds = inds(mod(1:size(inds,1),8)==5);
th_y_inds = inds(mod(1:size(inds,1),8)==6);
k_y_inds = inds(mod(1:size(inds,1),8)==7);
dR_y_inds = inds(mod(1:size(inds,1),8)==0);

%R0 = square
side_n = round(n/4);
side_len = (side_n+1)*ds;
R_0 = zeros(n,2);
for i = 1:side_n
    R_0(i,1) = i*ds;
    R_0(i+side_n,1) = side_len;
    R_0(i+2*side_n,1) = side_len - i*ds;
    R_0(i+3*side_n,1) = 0;
    R_0(i,2) = 0;
    R_0(i+side_n,2) = i*ds;
    R_0(i+2*side_n,2) = side_len;
    R_0(i+3*side_n,2) = side_len - i*ds;
end
[th_0,k_0] = getShape(R_0,ds);
size(Q_0(R_x_inds))
size(R_0(:,1))
Q_0 = zeros(8*n,1);
Q_0(R_x_inds) = R_0(:,1);
Q_0(R_y_inds) = R_0(:,2);
Q_0(th_x_inds) = th_0(:,1);
Q_0(th_y_inds) = th_0(:,2);
Q_0(k_x_inds) = k_0(:,1);
Q_0(k_y_inds) = k_0(:,2);

dQ_ode = @(t,y)dQ(t,y,ds,R_x_inds,th_x_inds,k_x_inds,dR_x_inds,R_y_inds,th_y_inds,k_y_inds,dR_y_inds);
M_ode = @(t,y)mass_mat(t,y,rho,@(k)k_t(k),ds,R_x_inds,th_x_inds,k_x_inds,dR_x_inds,R_y_inds,th_y_inds,k_y_inds,dR_y_inds);

[t,y] = ode45(dQ_ode,...
              [0 5],...
              Q_0,...
              odeset('Mass',M_ode));
          
plot3(y(:,1),y(:,2),t)
zlim([0 t(end)])
ylim([-2 12])
xlim([-2,12])

function ret = dQ(t,R,ds,R_x_inds,th_x_inds,k_x_inds,dR_x_inds,R_y_inds,th_y_inds,k_y_inds,dR_y_inds)
    t
    n = size(R_x_inds,1);
    R_x = R(R_x_inds);
    th_x = R(th_x_inds);
    k_x = R(k_x_inds);
    dR_x = R(dR_x_inds);
    R_y = R(R_y_inds);
    th_y = R(th_y_inds);
    k_y = R(k_y_inds);
    dr_y = R(dR_y_inds);
    
    A = zeros(size(R,1));
    %right half of binding velocities to derivative of positions
    i_minus = [n,1:n-1];
    i_plus = [2:n,1];
    for i = [i_minus;1:n;i_plus]
        %output row       %vector
       %bind derivative to pos
       A(R_x_inds(i(2)),dR_x_inds(i(2))) = 1;
       A(R_y_inds(i(2)),dR_y_inds(i(2))) = 1;
       %bind dRds to spacial derivative
       A(th_x_inds(i(2)),[th_x_inds(i(2)) R_x_inds(i(1)) R_x_inds(i(3))]) = [2,-1/ds,1/ds];
       A(th_y_inds(i(2)),[th_y_inds(i(2)),R_y_inds(i(1)),R_y_inds(i(3))]) = [2,-1/ds,1/ds];
       %bind d2Rds2 to second spacial derivative
       A(k_x_inds(i(2)),[k_x_inds(i(2)),R_x_inds(i(1)),R_x_inds(i(2)),R_x_inds(i(3))]) = [1,-1/(ds^2),2/(ds^2),-1/(ds^2)];
       A(k_y_inds(i(2)),[k_y_inds(i(2)),R_y_inds(i(1)),R_y_inds(i(2)),R_y_inds(i(3))]) = [1,-1/(ds^2),2/(ds^2),-1/(ds^2)];
       if i(2) > 12 && i(2) <= n-12 
            %A([R_x_inds(i(2))-12:R_x_inds(i(2))+12],[R_x_inds(i(2))-12:R_x_inds(i(2))+12])
            %i
       end
    end
    ret = A*R;
end

function M = mass_mat(t,R,rho,k_t,ds,R_x_inds,th_x_inds,k_x_inds,dR_x_inds,R_y_inds,th_y_inds,k_y_inds,dR_y_inds)
    t
    n = size(R_x_inds,1);
    R_x = R(R_x_inds);
    th_x = R(th_x_inds);
    k_x = R(k_x_inds);
    dR_x = R(dR_x_inds);
    R_y = R(R_y_inds);
    th_y = R(th_y_inds);
    k_y = R(k_y_inds);
    dR_y = R(dR_y_inds);
    
    M = zeros(size(R,1));
    for i = 1:n
        %left half of binding velocities to derivative of positions
        M(R_x_inds(i),R_x_inds(i)) = 1;
        M(R_y_inds(i),R_y_inds(i)) = 1;
        %EOM
        C = k_t(0)*cross([k_x(i),k_y(i),0],[th_x(i),th_y(i),0])*[0; 0; 1];
        M(dR_x_inds(i),[dR_x_inds(i),dR_y_inds(i),k_x_inds(i),k_y_inds(i),th_x_inds(i),th_y_inds(i)]) = [rho*dR_x(i),rho*dR_y(i),C*th_y(i),-C*th_x(i),-C*k_y(i),C*k_x(i)]; 
        if i > 12 && i <= n-12 
           M([R_x_inds(i)-12:R_x_inds(i)+12],[R_x_inds(i)-12:R_x_inds(i)+12])
           i
        end
    end
end

function ret = k_t(k)
    E = 10; %placeholders
    I = 5;  %placeholders
    ret = E*I;
end

function [th,k] = getShape(R,ds)
    %R is nx2 [R_x,R_y]
    n = size(R,1);
    th = zeros(size(R));
    k = zeros(size(R));
    i_minus = [n,1:n-1];
    i_plus = [2:n,1];
    for i = [i_minus;1:n;i_plus]
       th(i(2),1) = (R(i(3),1) - R(i(1),1))/(2*ds);
       th(i(2),2) = (R(i(3),1) - R(i(1),1))/(2*ds);
       k(i(2),1) = (R(i(1),1) - 2*R(i(2),1) + R(i(3),1))/(ds^2);
       k(i(2),2) = (R(i(1),2) - 2*R(i(2),2) + R(i(3),2))/(ds^2);
    end
end
