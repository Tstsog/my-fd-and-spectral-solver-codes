function [] = heat_equation_1D_btcs_scheme
%%% --- 1D heat equation with Dirichlet boundary condition
%
% du/dt = D*d^2u/dx^2 with Drichlet BC: u(a,t) = u(b,t) = 0;
%                                       u(x,0) = 2*sin(2*pi*x)
% analytic solution: u(x,t) = 2*exp(-(pi^2/4)*t)*sin(2*pi*x).
%
% The finite difference: the backward-time centered space (BTCS), that is,
%  matrix equation: 
% (I-lambda*A)*w^(n+1) = w^(n) + lambda*b^(n+1).
% A is tridiagonal matrix.
% (I-lambda*A) is evolution matrix.
% b^(n+1) = [u^(n+1)_{a} 0,...0, u^(n+1)_{b}]^(T), (T = transpose). 
% b^(n+1) = 0 for all n since u(a,t) = (b,t) = 0
% lambda = D*dt/dx^2.
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
N = 128;  % number of grid point of x axis
dx = (b-a)/N; % step size in x
% ---
x = zeros(N+1,1); % total number of points is N+1
for i = 1:N+1
    x(i) = a + (i-1)*dx;
end
% ---
%%% time parameter
ti = 0.; % t(0)
tf = 1.; % t(N)
Nt = 256; % number of grid point in time t
dt = (tf-ti)/Nt; % step size in time t
t = zeros(Nt,1); % total number of points 
%
lambda = D*dt/dx^2 ;
%  matrix equation is
% (I-lambda*A)*w^(n+1) = w^(n) + lambda*b^(n+1)
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
evolution_mat = unit_I - u_mat;
w0 = w0(2:N); % taking into account the boundary condition  
%
%%% time loop starts ---
for j = 1:Nt
    t(j) = ti + j*dt;
%
    w_old = w0;
    w_new = evolution_mat(2:N,2:N)\w_old; % % BTCS scheme: solves linear algebra problem (I-lambda*A)*w^(n+1) = w^(n) at each time step
    w0 = w_new;
%
end
%%% time loop ends ---
w_sol = [0.; w_new; 0.]; % numerical solution after taking into account the BC
%%% 
% ploting numerical and analytic solutions together
figure(1)
hold on
plot(x, w_sol, 'b')                                   % numerical solution
plot(x, 2.*exp(-(pi^2/4)*t(Nt)).*sin(2.*pi.*x), 'ro') % analytic solution
hold off
%axis([-25. 10. -1.100 0.20])
set(gca,'FontSize',18)
xlabel('x') % ,'fontsize',16
ylabel('u(x,1)') % ,'Rotation', 1
%%%

%return
end