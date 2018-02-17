function [cov_x,cov_y,cov_xy,Ahat,Bhat, Chat, Chat_star, Dhat,Ehat]=est_fun(n,h,grid_length,data)
% function [cov_x,cov_y,cov_xy,cov_y_mpp,expect_y,Ahat,Bhat, Chat, Chat_star, Dhat,Ehat]=est_fun(n,h,grid_length,data)
%cov_y_mpp: original mpp estimation for E(y(u)y(v))
%expect_y:   original mpp estimation for E(y(u))
%n=200; h=0.02;grid_length=100
T=1;
grid=T/grid_length/2:T/grid_length:T;
int_u=cdf('Normal',(1-grid)/h,0,1)-cdf('Normal',(-grid)/h,0,1);%1*50
edge=int_u'*int_u;


%  clear delta1 delta_t1 delta_t2 delta_t3 delta_t4
for i=1:n
%      i=1
    sub_data=data(data(:,1)==i,:);
    delta_sub=exp(-0.5*(bsxfun(@minus, sub_data(:,2),grid)/h).^2)/sqrt(2*pi)/h;%Ni*50
    delta(i,:)=sum(delta_sub); %Ni*50
    delta_sub_exp=kron(exp(sub_data(:,3)),ones(1,grid_length)).*delta_sub;
    
    delta1_sub_exp=kron((sub_data(:,3)),ones(1,grid_length)).*delta_sub;%for origianl marked PP estimation
    delta1_exp(i,:)=sum(delta1_sub_exp); %n*50;
    
    delta_exp(i,:)=sum(delta_sub_exp); %n*50;
    delta_c2(:,:,i)=delta_sub_exp'*delta_sub;
    delta_c2_star(:,:,i)=delta_sub'*delta_sub_exp;
    delta_d2(:,:,i)=delta_sub'*delta_sub;
    delta_e2(:,:,i)=delta_sub_exp'*delta_sub_exp;
    delta1_e2(:,:,i)=delta1_sub_exp'*delta1_sub_exp;
end
Ahat=sum(delta)./int_u/n;%1*50
Bhat=sum(delta_exp)./int_u/n;%1*50
C1hat=delta_exp'*delta;
C1hat_star=delta'*delta_exp;
C2hat=sum(delta_c2,3);
C2hat_star=sum(delta_c2_star,3);
D1hat=delta'*delta;
D2hat=sum(delta_d2,3);
E1hat=delta_exp'*delta_exp;
E2hat=sum(delta_e2,3);

Bhat_mpp=sum(delta1_exp)./int_u/n;%1*50
E1hat_mpp=delta1_exp'*delta1_exp;
E2hat_mpp=sum(delta1_e2,3);
Ehat_mpp=(E1hat_mpp-E2hat_mpp)./edge/n;

Chat=(C1hat-C2hat)./edge/n;
Chat_star=(C1hat_star-C2hat_star)./edge/n;
Dhat=(D1hat-D2hat)./edge/n;
Ehat=(E1hat-E2hat)./edge/n;


exp_covx=Dhat./kron(Ahat',ones(1,grid_length))./kron(ones(grid_length,1),Ahat);
exp_covxy=Chat.*kron(Ahat',ones(1,grid_length))./kron(Bhat',ones(1,grid_length))./Dhat;
exp_covy=Ehat.*Dhat./Chat./Chat_star ;
cov_x=log(exp_covx);cov_xy=log(exp_covxy);cov_y=log(exp_covy);
% cov_y_mpp=(Ehat_mpp./Dhat);
% expect_y=Bhat_mpp./Ahat;

