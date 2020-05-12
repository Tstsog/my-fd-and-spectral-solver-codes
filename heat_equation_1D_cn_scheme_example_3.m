function [] = heat_equation_1D_cn_scheme_example_3
%%% --- 1D heat equation 
%
% du/dt = D*d^2u/dx^2 - t + 8*x^2 with BC: u(0,t) = 0.
%                                          u(1,t) = 8*t
%                                          u(x,0) = 2*sin(2*pi*x)
% analytic solution: u(x,t) = 2*e^(-pi^2*t/4)sin(2*pi*x) + 8*x^2*t.
%
% The finite difference: the Crank-Nicolson scheme (CN), that is,
%  matrix equation: 
% (I-lambda*A)*w^(n+1) = (I+lambda*A)*w^(n) + lambda*(b^(n) + b^(n+1)) + 0.5*dt*(s^(n)+s^(n+1))
% A is tridiagonal matrix.
% (I-lambda*A) is evolution matrix.
% b^(n+1) = [u^(n+1)_{a} 0,...0, u^(n+1)_{b}]^(T), (T = transpose). 
% b^(n+1) for all n since u(a,t) = 0 & u(b,t) = 8*t
% lambda = D*dt/(2*dx^2).
% Reference book: D. Bradie, A friendly introduction to numerical analysis 
%
% Written by Tsogbayar Tsednee (PhD), California State University Northridge 
% Contact: tsog215@gmail.com
% Date: December 23, 2018
%%%
% ----------------------------------------------------------------
format short; clear; clc;
% ---
% grid and initial condition
D = 1./16.; % diffision parameter
a = 0.; % x(0)
b = 1.; % x(N+1)
N = 64;  % number of grid point of x axis
dx = (b-a)/N; % step size in x
% ---
x = zeros(N+1,1); % total number of points is N+1
for i = 1:N+1
    x(i) = a + (i-1)*dx;
end
% ---
%%% time parameter
ti = 0.; % t(0)
tf = 8.; % t(N)
Nt = 256; % number of grid point in time t
dt = (tf-ti)/Nt; % step size in time t
%
lambda = D*dt/(2.*dx^2) ;
%  matrix equation is
% (I-lambda*A)*w^(n+1) = (I+lambda*A)*w^(n) + lambda*(b^(n) + b^(n+1)) + 0.5*dt*(s^(n)+s^(n+1))
u_mat = zeros(N+1,N+1);
for i = 2:N
    u_mat(i,i-1) = lambda;
    u_mat(i,i) = -2.*lambda;
    u_mat(i,i+1) = lambda;
end
u_mat(1,1) = -2.*lambda; u_mat(N+1,N+1) = -2.*lambda ; 
u_mat(1,2) = lambda; u_mat(N+1,N) = lambda;
%u_mat;
%%%
unit_I = eye(N+1); % unit diagonal matrix 
% at t0 = 0
w0 = 2.*sin(2.*pi.*x); % % at t0 = 0.00 & initial
% evolution matrix 
evolution_mat_l = unit_I - u_mat;
evolution_mat_r = unit_I + u_mat;
w0 = w0(2:N); % taking into account the boundary condition  
%
tn = zeros(Nt,1); tnp1 = zeros(Nt,1); % time points
bn = zeros(Nt,1); bnp1 = zeros(Nt,1); % b(t(n)) & b(t(n+1)) vectors
%
%%% time loop starts ---
for j = 1:Nt
    tn(j) = ti + (j-1)*dt; % t(n)
    bn(j) = 8.*tn(j);    
    b_n = [0.; zeros(N-3,1); bn(j)];
    s_n = (-tn(j) + 8.*x(2:N).^2);
%
    tnp1(j) = ti + (j-0)*dt; % t(n + 1) 
    bnp1(j) = 8.*tnp1(j);    
    b_np1 = [0.; zeros(N-3,1); bnp1(j)];
    s_np1 = (-tnp1(j) + 8.*x(2:N).^2);
%
    w_old = w0;
    w_r_new = evolution_mat_r(2:N,2:N)*w_old + lambda.*(b_n + 1.*b_np1) + 0.5*dt.*(s_n+s_np1); % 
    w_new = evolution_mat_l(2:N,2:N)\w_r_new; % 
    w0 = w_new ;    
%
end
%%% time loop ends ---
time = tnp1(Nt);
%%%
w_num = [0.; w_new; 8.*time ];
u_exact = 2.*exp(-pi^2*time/4).*sin(2*pi*x) + 8.*x.^2.*time;

% ploting numerical and analytic solutions together
figure(1)
hold on
plot(x, w_num, 'b')       % plot numerical solution 
plot(x, u_exact, 'ro')    % plot analytical solution
hold off
set(gca,'FontSize',18)
xlabel('x') % ,'fontsize',16
ylabel('u(x,1)') % ,'Rotation', 1
%
figure(2)
plot(x, w_num-u_exact) % plot difference between numerical and analytical solution
set(gca,'FontSize',16)
xlabel('x') % ,'fontsize',16
ylabel('u_{num}-u_{exact}') % ,'Rotation', 1

end
