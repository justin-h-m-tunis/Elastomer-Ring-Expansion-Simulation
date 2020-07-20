function [th,k,j,dsds] = getShape(R,ds,use_scale)
    if nargin == 2
        use_scale = true;
    end
    %R is nx2 [R_x,R_y]
    n = size(R,1);
    th = zeros(size(R));
    k = zeros(size(R));
    j = zeros(size(R));
    dsds = zeros(size(R));
    i_mminus = [n-1,n,1:n-2];
    i_minus = [n,1:n-1];
    i_plus = [2:n,1];
    i_pplus = [3:n,1,2];
    for i = [i_mminus;i_minus;1:n;i_plus;i_pplus]
       dth = (-1*R(i(5),:) + 8*R(i(4),:) - 8*R(i(2),:) + R(i(1),:))/12;
       dth_b = (R(i(3),:) - R(i(1),:))/2;
       dth_f = (-R(i(3),:) + R(i(5),:))/2;
       dk = (-1*R(i(1),:) + 16*R(i(2),:) - 30*R(i(3),:) + 16*R(i(4),:) - R(i(5),:))/12;
       DS = norm(dth);%adjusted ds
       DS_b = norm(dth_b);
       DS_f = norm(dth_f);
       th(i(3),:) = dth./DS;
       k(i(3),:) = dk./(DS^2);
       dsds(i(3),:) = dth_b.*((DS_b-ds)/DS_b) - dth_f.*((DS_f-ds)/DS_f);
       j(i(3),1) = (-R(i(1),1) + 2*R(i(2),1) - 2*R(i(4),1) + R(i(5),1))/(2*DS^3);
       j(i(3),2) = (-R(i(1),2) + 2*R(i(2),2) - 2*R(i(4),2) + R(i(5),2))/(2*DS^3);
    end
end