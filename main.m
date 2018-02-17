clear all
%% simulation
T=1;
n=200;
rep=100;
eigen_values=[0.5,0.2];
lambda0=100;
mu0=0;
parfor i=1:rep
    %     i
    %x, y are not independent
    mu=zeros(1,4);
    sigma1=[0.5 0; 0 0.2];
    sigma2=[0.2 0.15;0.15 0.1]; 
    %     sigma2=[0 0; 0 0];% if x, y are independent
    sigma_m=[sigma1 sigma2;sigma2' sigma1];
    
    xi=mvnrnd(mu,sigma_m,n);
    [data]=sim_fun(lambda0,mu0,n,T,xi);
    DATA{i}=data;
    XI(:,:,i)=xi;
end
save mfpp_n200_dependent_error DATA XI mu0 lambda0
%% estiamtion
T=1;grid_length=100;
grid=T/grid_length/2:T/grid_length:T;
h=0.02;n=200;rep=100;
tic
parfor i=1:rep
     i
    [cov_x(:,:,i),cov_y(:,:,i),cov_xy(:,:,i),Ahat(:,:,i),Bhat(:,:,i), Chat(:,:,i), ...
        Chat_star(:,:,i), Dhat(:,:,i),Ehat(:,:,i)]=est_fun(n,h,grid_length,DATA{i});
end
toc
% save est_quantities_h002 cov_x cov_y cov_xy Ahat Bhat Chat Chat_star Dhat Ehat

%eigen-decompostion for many simulations
% load est_quantities_h002
parfor i=1:rep
    % traditional way to estiamte marked point process
    %     expect_yy(:,:,i)=expect_y(:,:,i)'*expect_y(:,:,i);
    %     cov_xy_mpp(:,:,i)=cov_y_mpp(:,:,i)-expect_y(:,:,i)'*expect_y(:,:,i);
         
    [V_x,D_x]=eig((cov_x(:,:,i)+cov_x(:,:,i)')/2/grid_length);
    [V_y,D_y]=eig((cov_y(:,:,i)+cov_y(:,:,i)')/2/grid_length);
    [V_xy,D_xy]=eig((cov_xy(:,:,i)+cov_xy(:,:,i)')/2/grid_length);
    Dx_all(:,i)=flip(diag(D_x));   Dy_all(:,i)=flip(diag(D_y));  Dxy_all(:,i)=flip(diag(D_xy));
    Vx_all(:,:,i)=flip(grid_length^0.5*V_x,2);
    Vy_all(:,:,i)=flip(grid_length^0.5*V_y,2);
    Vxy_all(:,:,i)=flip(grid_length^0.5*V_xy,2);
end
mean_all=[mean(squeeze(Dx_all(1,:,:))),    mean(squeeze(Dx_all(2,:,:))),  ...
    mean(squeeze(Dy_all(1,:,:))),    mean(squeeze(Dy_all(2,:,:))),...
    mean(squeeze(Dxy_all(1,:,:))),    mean(squeeze(Dxy_all(2,:,:)))];

var_all=[std(squeeze(Dx_all(1,:,:))),std(squeeze(Dx_all(2,:,:))), std(squeeze(Dy_all(1,:,:))),...
    std(squeeze(Dy_all(2,:,:))),     std(squeeze(Dxy_all(1,:,:))),std(squeeze(Dxy_all(2,:,:)))];
% [mean_all var_all]

eigenvalue=table({'X1';'X2';'Y1';'Y2';'XY1';'XY2'},[mean_all'], [var_all']);
eigenvalue.Properties.VariableNames = {'direction' 'mean' 'std'};
eigenvalue

figure
subplot(2,2,1);plot(squeeze(Vx_all(:,1,:)));ylim([-2,2])
subplot(2,2,2);plot(squeeze(Vx_all(:,2,:)));ylim([-2,2])
subplot(2,2,3);plot(squeeze(Vy_all(:,1,:)));ylim([-2,2])
subplot(2,2,4);plot(squeeze(Vy_all(:,2,:)));ylim([-2,2])
