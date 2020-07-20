function [th, k,dsds,k_dot, dsds_dot] = get_dshapedt(R,dR,ds)
    %R is nx2 [R_x,R_y]
    n = size(R,1);
    th = zeros(size(R));
    k = zeros(size(R));
    dsds = zeros(size(R));
    k_dot = zeros(size(dR));
    dsds_dot = zeros(size(dR));
    i_mminus = [n-1,n,1:n-2];
    i_minus = [n,1:n-1];
    i_plus = [2:n,1];
    i_pplus = [3:n,1,2];
    for i = [i_mminus;i_minus;1:n;i_plus;i_pplus]
       dth = (-1*R(i(5),:) + 8*R(i(4),:) - 8*R(i(2),:) + R(i(1),:))/12;
       dth_dot = (-1*dR(i(5),:) + 8*dR(i(4),:) - 8*dR(i(2),:) + dR(i(1),:))/12;
       dth_b = (R(i(3),:) - R(i(1),:))/2;
       dth_b_dot = (dR(i(3),:) - dR(i(1),:))/2;
       dth_f = (-R(i(3),:) + R(i(5),:))/2;
       dth_f_dot = (-dR(i(3),:) + dR(i(5),:))/2;
       dk = (-1*R(i(1),:) + 16*R(i(2),:) - 30*R(i(3),:) + 16*R(i(4),:) - R(i(5),:))/12;
       dk_dot = (-1*dR(i(1),:) + 16*dR(i(2),:) - 30*dR(i(3),:) + 16*dR(i(4),:) - dR(i(5),:))/12;
       DS = norm(dth);%adjusted ds
       DS_dot = (dth*dth_dot.')/DS;
       th(i(3),:) = dth./DS;
       k(i(3),:) = dk./(DS^2);
       k_dot(i(3),:) = (dk_dot.*DS - dk.*DS_dot)./(DS^2);
       DS_b = norm(dth_b);
       DS_b_dot = (dth_b*dth_b_dot.')/DS_b;
       DS_f = norm(dth_f);
       DS_f_dot = (dth_f*dth_f_dot.')/DS_f;
       dsds(i(3),:) = dth_b.*(1-ds/DS_b) - dth_f.*(1-ds/DS_f);
       dsds_dot(i(3),:) = dth_b_dot.*(1-ds/DS_b) + dth_b.*(ds*DS_b_dot*DS_b^(-2)) - dth_f_dot.*(1-ds/DS_f) + dth_f.*(ds*DS_f_dot*DS_f^(-2));
    end
end