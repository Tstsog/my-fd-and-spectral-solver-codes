%%%
%%% --- 1D stochastic differential equation - Ornstein-Uhlenbeck process
% dx(t) = lambda*x(t)*dt + mu*x(t)*dw(r) & x(0) = x0 & 0 < t < T
% analytic solution: x(t) = x0*exp((lambda - mu^2/2).*t + mu*w(t)),
%
% The finite difference is Euler's method. 
%%% ---
%
% Reference: SIAM Review v43, p526 (2001) 
%
% Written by Tsogbayar Tsednee (PhD), Kyungpook National University, South Korea 
% Contact: tsog215@gmail.com
% Date: Jan 11, 2021
%
%
function [] = stochastic_diff_eq_euler_method_example_1
%
clear;
clc;
randn('state', 100)
%
N = 2000.;                      % number of the grid point 
T = 10.;                        % T is end time 
dt = T/N;                       % dt is time step
x0 = 10.0;                      % initial value x(0) at t=0
lambda = 0.50; 
mu = 0.80;
%
t = [dt:dt:T];
%
dw = sqrt(dt)*randn(1,N);        %  Brownian increments
w = cumsum(dw,2);                %  Discritized Brownian path
%
x_exact = x0.*exp((lambda - mu^2/2).*t + mu.*w);  % analytical solution
%
x = zeros(N,1);
x(1) = x0;
for i = 1:N
   x(i+1) = x(i) + lambda*dt*x(i) + mu*x(i)*dw(i);  % Euler's scheme
end
x';
figure(1)
hold on
plot([0:dt:T], [x0,x_exact], 'ro','LineWidth', 1.2) % exact
plot([0:dt:T], x, 'b-','LineWidth', 1.2)
hold off
%axis([0 1 0 1.2])
xlabel('t')
ylabel('x(t)', 'Rotation', 90, 'VerticalAlignment','middle')
legend('exact', 'numer')
set(gca,'FontSize',18)









