%%%
%%% --- 1D stochastic differential equation - Ornstein-Uhlenbeck process
% dx(t) = -x(t)*dt/(2*tau) + dw(r)/sqrt(tau) & x(0) = x0 & 0 < t < T
% analytic solution: x(t) = x0*exp(-t/(2*tau)) + sqrt(1-exp(-t/tau))*w(t),
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
function [] = stochastic_diff_eq_euler_method_example_2
%
clear all;
clc;
randn('state', 100)
%
N = 1000.;                         % number of grid points
T = 1.0;                           % end time
dt = T/N;                          % time step
x0 = 0.5;                          % initial value
tau = 1.5;
t = [dt:dt:T];
%
dw = sqrt(dt)*randn(1,N);          %  Brownian increments
w = cumsum(dw,2);                  %  Discritized Brownian path
%
x_exact = x0.*exp(-t./(2.*tau)) + 1.*sqrt(1.-exp(-t./tau)).*w; % analytical solution
%
x = zeros(N,1);
x(1) = x0;
for i = 1:N
   x(i+1) = x(i) - dt*x(i)/(2*tau) + dw(i)/sqrt(tau);   % Euler's scheme
end
x';
figure(1)
hold on
plot([0:dt:T], [x0,x_exact], 'r','LineWidth', 1.2) % exact
plot([0:dt:T], x, 'b-','LineWidth', 1.2)
hold off
%axis([0 1 0 1.2])
xlabel('t')
ylabel('x(t)', 'Rotation', 90, 'VerticalAlignment','middle')
legend('exact', 'numer')
set(gca,'FontSize',18)

