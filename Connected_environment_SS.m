clc
clear all
close all
fprintf('Code is running...')
%% Loading input data and essential files
load('connSS_60000_trj1.mat') %Trajectory
load('connSS_60000_speed_trj1.mat') %Speed
load('connSS_60000_hd_trj1.mat') %Head direction
load('hd_s1_s2_wt_som') %SOM weights of Head Direction layer
%% Plotting the Trajectory
figure; plot(pos(:,1),pos(:,2))
axis off
title('Trajectory of Square-Square connected environment')
%% Head Direction (HD) response computing
speed = speed;
phase1d = zeros(100,1);
PI2d = zeros(10);
X=[];
trj_hd_resp=[];
for ii = 1:size(speed,1)
    X1 = [cosd(theta_real_deg(1)) sind(theta_real_deg(1))]; X2 = [cosd(theta_real_deg(ii)) sind(theta_real_deg(ii))];
    s1 = X2(1)*X1(2) - X1(1)*X2(2);  %sin(theta1-theta2)
    s2 = X2(1)*X1(1) + X1(2)*X2(2);  %cos(theta1-theta2)
    X=[s1 s2];
    y = respsom2dlinear(X,wt2);
    trj_hd_resp(:,:,ii) = y;    
end
%% Path Integration (PI) oscillations
X = zeros(100,1); Y = ones(100,1); %Initializing the State variables
dt = 0.01;
bf = 6*2*pi; %Base Frequency of the PI oscillations
niter = size(trj_hd_resp,3);
beta = 50; %Spatial scaling parameter
t = 0; 
Xbg = 1; Ybg = 0;
tarr=[];
theta = zeros(100,1);
for ii = 2:niter
    y = trj_hd_resp(:,:,ii); 
    inp1d = reshape(y,100,1);
    thetadot = bf + beta*speed(ii)*inp1d*10;
    theta(:,ii)=theta(:,ii-1)+thetadot*dt;
end
Xarr=cos(theta);
PI1d=Xarr; 
%% Lateral Anti-Hebbian network (LAHN)
[N K] = size(PI1d); %N --> Dimension    K---> # of samples
PI1d=removemean(PI1d);
alphaa = 0.000001/K; %Learning rate of the afferent connections  
betaa = 0.000001/K; %Learning rate of the lateral connections 
output_neuron_nmbr = 13; %No.of neurons in the LAHN
maxiter = 2000000; %Maximum training iterations
%% Single cell firing field during the initial phase of LAHN training
load('T_initial')  %LAHN weights collected during the initial phase of training
Xarr = cos(theta);
PI1d = Xarr;
neuron_num = 8; %Select the neuron number
w=T(neuron_num,:);w = w'; %Weights from PI to the selected LAHN neuron
ot=w'*PI1d; ot=ot'; %Output activity of the LAHN neuron
thresh=max(abs(ot))*.65;
firr=find(abs(ot)>thresh);
firposgrid=pos(firr,:); %Firing positions on the trajectory
figure; plot(pos(:,1),pos(:,2)); hold on; plot(firposgrid(:,1),firposgrid(:,2),'.r', 'markersize', 15); 
axis off
title('Firing field of the neuron during the initial phase of training')
%% Single cell firing field during the mid phase of LAHN training
load('T_mid') %LAHN weights collected during the mid phase of training
Xarr = cos(theta);
PI1d = Xarr;
neuron_num = 8;
w=T(neuron_num,:);w = w';
ot=w'*PI1d; ot=ot';
thresh=max(abs(ot))*.65;
firr=find(abs(ot)>thresh);
firposgrid=pos(firr,:);
figure; plot(pos(:,1),pos(:,2)); hold on; plot(firposgrid(:,1),firposgrid(:,2),'.r', 'markersize', 15); 
axis off
title('Firing field of the neuron during the mid phase of training')
%% Single cell firing field post LAHN training
load('T_final') %LAHN weights collected post training
Xarr = cos(theta);
PI1d = Xarr;
neuron_num = 8;
w=T(neuron_num,:);w = w';
ot=w'*PI1d; ot=ot';
thresh=max(abs(ot))*.65;
firr=find(abs(ot)>thresh);
firposgrid=pos(firr,:);
figure; plot(pos(:,1),pos(:,2)); hold on; plot(firposgrid(:,1),firposgrid(:,2),'.r', 'markersize', 15); 
axis off
title('Firing field of the neuron post training')
fprintf('Simulation complete')