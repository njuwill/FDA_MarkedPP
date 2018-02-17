% this is the function to simulate marked temporal point processes
%eigen function for X is (1 sqrt(2)*sin(2*pi*t))
%eigen functions for Y are 1 and sqrt(2)*sin(4*pi*t)

%


function [data]=sim_fun(lambda0,mu0,n,T,xi)

xi_x_1=xi(:,1);xi_x_2=xi(:,2);xi_y_1=xi(:,3);xi_y_2=xi(:,4);

%now simulate the events according to the maximum intensities for each
%subject and on each day

% max_intensity=lambda0*exp(sqrt(xi_x_1.^2+xi_x_2.^2)*sqrt(2));  %one setting
max_intensity=lambda0*exp(xi_x_1+abs(xi_x_2)*sqrt(2));         %another setting
counts_before_thinning=random('Poisson',max_intensity);        %number of events to be simulated according the maximum intensity


counts_before_thinning_nonzero_id=find(counts_before_thinning>0);
counts_before_thinning_nonzero=counts_before_thinning(counts_before_thinning_nonzero_id);
counts_before_thinning_total=sum(sum(counts_before_thinning_nonzero));

%now find out the subjects and days for each nonzero count
subjects_all=[1:n];
subjects_nonzero=subjects_all(counts_before_thinning_nonzero_id);
max_intensity_nonzero=max_intensity(counts_before_thinning_nonzero_id);

starting_id=1;
subjects_nonzero_new=[];max_intensity_nonzero_new=[];
for i=1:length(counts_before_thinning_nonzero_id)
    ending_id=starting_id+counts_before_thinning_nonzero(i)-1;
    subjects_nonzero_new(starting_id:ending_id,1)=subjects_nonzero(i);
    max_intensity_nonzero_new(starting_id:ending_id,1)=max_intensity_nonzero(i);
    starting_id=ending_id+1;
end

events_before_thinning=random('Uniform',0,T,counts_before_thinning_total,1);
prob_thinning=random('Uniform',0,1,counts_before_thinning_total,1);

intensity=lambda0.*exp(xi_x_1(subjects_nonzero_new)+xi_x_2(subjects_nonzero_new).*sin(2*pi*events_before_thinning/(T))*sqrt(2));
% another setting
% intensity=lambda0.*exp(xi_x_1(subjects_nonzero_new).*sin(2*pi*events_before_thinning)*sqrt(2)+...
%     xi_x_2(subjects_nonzero_new).*cos(2*pi*events_before_thinning)*sqrt(2));


id_kept=find(intensity./max_intensity(subjects_nonzero_new)>=prob_thinning);
t=events_before_thinning(id_kept); %N*1
subject_id=subjects_nonzero_new(id_kept);
%first way to get zt
% yt=cdf('Normal',t,zeros(1,length(t)),sqrt(2)*abs(sin(4*pi*t))+0.1);
% yt_fun=@(t) cdf('Normal',t,0,sqrt(2)*abs(sin(4*pi*t))+0.1);
% zt= sqrt(2)*cos(4*pi*t)+yt_fun(t) ;

%second way to get zt
% zt=sqrt(2)*cos(4*pi*t)+xi_y_1(subject_id)+sqrt(2)*xi_y_2(subject_id).*sin(4*pi*t);
zt=mu0+xi_y_1(subject_id)+sqrt(2)*xi_y_2(subject_id).*sin(4*pi*t);
% zt=mu0+xi_y_1(subject_id)+sqrt(2)*xi_y_2(subject_id).*sin(4*pi*t)+random('Normal',0,0.2,length(t),1);
data=[subject_id,t,zt];
