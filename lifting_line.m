%% This Script uses the lifting line method to calculate induced drag

% Also caluculates yaw and roll moment, as well as the induced drag using

%% Purpose

% This scripts allows the lift distribution to be altered to fit different
% flight situations at the induced yaw and roll moment can be determined
% from it.

%% Assumptions

% Approximation based on Ludwing Prandtl's lifting line method
%(Kelvin and Helmholtz Vortex theorems combined with the Kutta-Zhukovsky
%theory of lift: Produces finite-wing affects)

%Inviscid Flow
%single lifting vortex

%% Clear workspace
clc
clear

%% Define wing geometry in meters
Uinf = 1;
b = 2;
max_chord = 0.3;
min_chord = 0.1;
mass = 1;
rho = 1;
A = -1/3;

% Exponent of the funtion 1/2 is an ellipse, 3/2 is the bell curve
exponent = 3/2;

%% Evaluate Properties

%number of points: choose an odd number
points = 10001;

dy = 1/(points - 1);

%Create the bell curve lift distribution
y = b * (-1/2:dy:1/2);
theta = acos(-2/b * y);

len = length(y);
Gamma = sin(theta) + A .* sin(3 * theta);
Gamma = Gamma/(max(Gamma));
% Use the forward-center-backwards difference method to calculate the 
% change in lift and change in drag
dGammady = FCB(Gamma,y,len);

% Calculate the induced angle of attack alpha
const = 1/(4 * pi * Uinf); %This will be multiplied by the integral
alpha = zeros(1,len);
for iindex = 1:len
    dalpha = zeros(1,len);
    for jindex = 1:len
        if jindex ~= iindex
            dalpha(jindex) = dGammady(jindex)/(y(iindex) - y(jindex));
        elseif iindex == 1
            dalpha(jindex) = (dGammady(jindex)/...
                (y(iindex) - y(jindex + 1)) + dGammady(jindex)/...
                (y(iindex + 1) - y(jindex)))/2;
        elseif iindex == len
            dalpha(jindex) = (dGammady(jindex)/...
                (y(iindex - 1) - y(jindex)) + dGammady(jindex)/...
                (y(iindex) - y(jindex - 1)))/2;
        else
            dalpha(jindex) = (dGammady(jindex)/...
                (y(iindex - 1) - y(jindex)) + dGammady(jindex)/...
                (y(iindex) - y(jindex - 1)))/2;
        end
    end
    %Simpson 1/3 rule
    for jindex = 2:2:len
        dalpha(jindex) = 4 * dalpha(jindex);
    end
    for jindex = 3:2:len - 2
        dalpha(jindex) = 2 * dalpha(jindex);
    end    
    alpha(iindex) = -const * 1/3 * dy * sum(dalpha);
end  

% Calculate the induced drag
induced_drag = sin(alpha) .* rho * Uinf .* Gamma;
drag = round(dy * sum(sin(alpha) .* rho * Uinf .* Gamma)/0.158443596484950,15)

%Calculate the roll moment
roll_moment = round(- dy * sum(y .* rho * Uinf .* Gamma),15);
yaw_moment = round(-dy * sum(induced_drag .* y),15);

%% Plot the data
%Lift Distribution
figure(1)
subplot(3,1,1)
plot(y,rho * Uinf * Gamma)
grid on
title('Bell Curved Lift Distribution')
ylabel('L(y) [dimless]')
xlabel('y [dimless]')
subplot(3,1,2)
plot(y,alpha)
grid on
ylabel('DW [dimless]')
xlabel('y [dimless]')
subplot(3,1,3)
plot(y,induced_drag)
ylabel('D_{i}(y)[dimless]')
xlabel('y [dimless]')
grid on

%% This function uses the forward center backwards difference method
function [dydx] = FCB(y,x,len)
dydx = zeros(1,len);
for index = 1:len
    if index == 1
        dydx(index) = (y(index + 1) - y(index))...
            /(x(index + 1) - x(index));
    elseif index == len
        dydx(index) = (y(index) - y(index - 1))...
            /(x(index) - (x(index-1)));
    else
        dydx(index) = (y(index + 1) - y(index - 1))...
            /(x(index + 1) - (x(index -1)));
    end
end
end
