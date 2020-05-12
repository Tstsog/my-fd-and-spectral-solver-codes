function [] = heat_equation_1D_cn_scheme_example_2
%%% --- 1D heat equation 
%
% du/dt = D*d^2u/dx^2 with Drichlet BC: u(a,t) = u(b,t) = 0.5*(1+e^(-4*t));
%                                       u(x,0) = cos(x)^2
% analytic solution: u(x,t) = 0.5+0.5*e^(-4*t)cos(2*x),
%
% The finite difference: the Crank-Nicolson scheme (CN), that is,
%  matrix equation: 
% (I-lambda*A)*w^(n+1) = (I+lambda*A)*w^(n) + lambda*(b^(n) + b^(n+1))
% A is tridiagonal matrix.
% (I-lambda*A) is evolution matrix.
% b^(n+1) = [u^(n+1)_{a} 0,...0, u^(n+1)_{b}]^(T), (T = transpose). 
% b^(n+1) for all n since u(a,t) = (b,t) = 0.5*(1+e^(-4*t))
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
D = 1.; % diffision parameter
a = 0.; % x(0)
b = pi; % x(N+1)
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
tf = 4.; % t(N)
Nt = 256; % number of grid point in time t
dt = (tf-ti)/Nt; % step size in time t
%
lambda = D*dt/(2.*dx^2) ;
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
    bn(j) = 0.5*(1.+exp(-4.*tn(j)));    
    b_n = [bn(j); zeros(N-3,1); bn(j)];
%
    tnp1(j) = ti + (j-0)*dt; % t(n + 1) 
    bnp1(j) = 0.5*(1.+exp(-4.*tnp1(j)));    
    b_np1 = [bnp1(j); zeros(N-3,1); bnp1(j)];
%
    w_old = w0;
    w_r_new = evolution_mat_r(2:N,2:N)*w_old + lambda.*(b_n + 1.*b_np1); % 
    w_new = evolution_mat_l(2:N,2:N)\w_r_new; % 
    w0 = w_new ;    
%
end
%%% time loop ends ---
time = tnp1(Nt);
%%%
w_sol = [0.5*(1.+exp(-4.*time)); w_new; 0.5*(1.+exp(-4.*time)) ];
%
% ploting numerical and analytic solutions together
figure(1)
hold on
plot(x, w_sol, 'b') % plot numerical solution
hold off
%axis([a b 0 0.50])
set(gca,'FontSize',18)
xlabel('x') % ,'fontsize',16
ylabel('u(x,1)') % ,'Rotation', 1

%%%   
u_exact = 0.5 + 0.5.*exp(-4.*time).*cos(2.*x);
%
figure(2)
plot(x, w_sol-u_exact, 'b') % plot difference between numerical and analytical solution
set(gca,'FontSize',16)
xlabel('x') % ,'fontsize',16
ylabel('u_{num}-u_{exact}') % ,'Rotation', 1

end

