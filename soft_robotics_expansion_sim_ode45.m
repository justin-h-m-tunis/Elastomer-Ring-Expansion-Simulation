%squiggly
R0 = .03
r0 = .01
q = 5
R = @(th,R0,r0,q)[(R0 + r0.*cos(q.*th)).*cos(th),(R0 + r0.*cos(q.*th)).*sin(th)]
dR = @(th,R0,r0,q)[(-q*r0*sin(q.*th))*cos(th) - (R0 + r0.*cos(q.*th))*sin(th),(-q*r0*sin(q.*th))*sin(th) + (R0 + r0.*cos(q.*th))*cos(th)]
S = integral(@(th)norm(dR(th,R0,r0,q)),0,2*pi,'ArrayValued',true)
n = 100
ds = S/n
thetas = zeros(n,1);
R_0 = zeros(n,2);
R_0(1,:) = R(0,R0,r0,q);
for i = 2:n
    thetas(i) = thetas(i-1) + fminbnd(@(th)abs(integral(@(theta)norm(dR(theta,R0,r0,q)),thetas(i-1),thetas(i-1) + th,'ArrayValued',true)-ds),0,2*pi,optimset('TolFun',1e-10,'TolX',1e-10));
    R_0(i,:) = R(thetas(i),R0,r0,q);
end
R_0
thetas
rho = 995*ds*.01*.002;
tan_del = .125
wn_a = sqrt(k_a(0)/rho)
wn_t = sqrt(k_t(0)/rho)
c_a = wn_a*tan_del
c_t = wn_t*tan_del

r_int_max = .05;

inds = [1:n*4].';
R_x_inds = inds(mod(1:size(inds,1),4)==1);
dR_x_inds = inds(mod(1:size(inds,1),4)==2);
R_y_inds = inds(mod(1:size(inds,1),4)==3);
dR_y_inds = inds(mod(1:size(inds,1),4)==0);
Q_0 = zeros(4*n,1);
Q_0(R_x_inds) = R_0(:,1);
Q_0(R_y_inds) = R_0(:,2);

P = @(t)P_ext(t,0,.25,3,.05,r_int_max,n)

dQ_ode = @(t,y)dQ(t,y,rho,c_t,c_a,@(k)k_t(k),@(d)k_a(d),ds,R_x_inds,dR_x_inds,R_y_inds,dR_y_inds,P);
J_pattern = repmat(spdiags(ones(2*n,10),-4:5,2*n,2*n),2,2)
%M_ode = @(t,y)mass_mat(t,y,rho,@(k)k_t(k),ds,R_x_inds,dR_x_inds,R_y_inds,dR_y_inds);
tspan = [0 4]

[t,y] = ode15s(dQ_ode, tspan,Q_0);


figure(1)
hold on
for i = 1:size(t,1)
    %plot3(y(i,R_x_inds),y(i,R_y_inds),t(i)*ones(size(R_x_inds,1),1),'blue')
    %[th,k,j,del_ds] = getShape([y(i,R_x_inds).',y(i,R_y_inds).'],ds);
    %quiver3(y(i,R_x_inds).',y(i,R_y_inds).',t(i)*ones(size(R_x_inds,1),1),k(:,1),k(:,2),zeros(size(R_x_inds,1),1),'red')
    %quiver3(y(i,R_x_inds).',y(i,R_y_inds).',t(i)*ones(size(R_x_inds,1),1),k(:,1),k(:,2),zeros(size(R_x_inds,1),1),'magenta')
end
zlim([0 t(end)])
ylim([-2*R0, 2*R0])
xlim([-2*R0, 2*R0])
hold off
while true
    expansion_animation(t,y,ds,rho,c_t,c_a,R_x_inds,R_y_inds,dR_x_inds,dR_y_inds,P,R0,2,18,tspan,'expansion.gif')
    expansion_animation(t,y,ds,rho,c_t,c_a,R_x_inds,R_y_inds,dR_x_inds,dR_y_inds,P,R0,2,50,tspan,'expansion.gif')
end
function ret = dQ(t,R,rho,c_t,c_a,k_t,k_a,ds,R_x_inds,dR_x_inds,R_y_inds,dR_y_inds,P_ext) %F_t should return 100x1 array of pressure
    t
    n = size(R_x_inds,1);
    R_x = R(R_x_inds);
    dR_x = R(dR_x_inds);
    R_y = R(R_y_inds);
    dR_y = R(dR_y_inds);
    [th,k,del_ds,k_dot,del_ds_dot] = get_dshapedt([R_x,R_y],[dR_x,dR_y],ds);
    P = P_ext(t);
    AQ = zeros(size(R));
    AQ(R_x_inds) = dR_x;
    AQ(R_y_inds) = dR_y;
    %include in report: ode settings
                       %fbd's
                       %explain c is viscoelastic, not surrounding fluid
                       %go into depth for delta s
    AQ(dR_x_inds) = k_t(0)/rho*k(:,1) + c_t*k_dot(:,1) - c_a*del_ds_dot(:,1) - k_a(0)/(rho)*del_ds(:,1) - P.'*ds.*R_x.*abs(th(:,2));
    AQ(dR_y_inds) = k_t(0)/rho*k(:,2) + c_t*k_dot(:,2) - c_a*del_ds_dot(:,2) - k_a(0)/(rho)*del_ds(:,2) - P.'*ds.*R_y.*abs(th(:,1));
    
    ret = AQ;
end

function P = P_ext(t,P_max,f_undulation,f_spacial,f_rotation,r_int_max,n)
    %rotating, undulating sin wave
    %P = (sin(f_undulation*t)+1)*P_max/2*(sin(2*pi*f_spacial/n*x + 2*pi*f_rotation) + 1);
    P = (cos(2*pi*f_undulation*t+pi)+1)*P_max/(2*r_int_max);
end
