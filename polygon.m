function[x,y]=polygon(num)
% close all
% clear all
% clc
% num =6; % Specify polygon variables
length = 1;
x = zeros(num,1); % preallocate x and y
y = zeros(num,1);
circleAng = 2*pi;% angle of the unit circle in radians
angleSeparation = circleAng/num; % the average angular separation between points in a unit circle
angleMatrix = 0: angleSeparation: circleAng; %create the matrix of angles for equal separation of points
angleMatrix(end) = []; % drop the final angle since 2Pi = 0

% generate the points x and y
for k = 1:num
    x1(k) = length * cos(angleMatrix(k));
    y1(k) = length * sin(angleMatrix(k));
end
x2=length * cos(angleMatrix(1));
y2=length * sin(angleMatrix(1));
x=[x1 x2];
y=[y1 y2];
%  figure; plot(x,y)
end






