


clc
clear all
close all
fprintf('Code is running...')
%% Loading input data and essential files
load('RECT_60000_hd_trj1.mat') %Trajectory
load('RECT_60000_speed_trj1.mat') %Speed
load('RECT_60000_trj1.mat') %Head direction
load('hd_s1_s2_wt_som') %SOM weights of Head Direction layer
%% Plotting the Trajectory
figure; plot(pos(:,1),pos(:,2))
axis off
title('Trajectory of the Transforming Environment')
