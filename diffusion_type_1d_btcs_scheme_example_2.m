function [] = diffusion_type_1d_btcs_scheme_example_2
%%%
%%% --- 1D convection-diffusion equation with Dirichlet boundary condition
% du/dt = -alpha*du/dx - 2*u + D*d^2u/dx^2 + f(x,t),
%                                            f(x,t) = (x^2 + 2*x - 2)*exp(-t)
%                                            u(0,t) = 0.0
%                                            u(1,t) = exp(-t)
%                                            u(x,0) = x^2
% analytic solution: u(x,t) = x^2*exp(-t),
% The finite difference we use here is the backward-time centered space (BTCS). 
%%% ---
% matrix equation: 
% (I + lambda1*u_x + dt*diag(u) - lambda2*u_xx)*w^(n+1) = w^(n) + 
%                                                         lambda_1*b1^(n+1) + 
%                                                         lambda_2*b2^(n+1)
%                                                         dt*s^(n+1) 
% A is tridiagonal matrix.
% (I + lambda1*u_x + dt*diag(u) - lambda2*u_xx) is evolution matrix.
% b1^(n+1) = [u^(n+1)_{a} 0,...0, -u^(n+1)_{b}]^(T), (T = transpose). 
% b1^(n+1) = 0 for all n since u(a,t) = (b,t) = 0
% b2^(n+1) = [u^(n+1)_{a} 0,...0, u^(n+1)_{b}]^(T), (T = transpose). 
% b2^(n+1) = 0 for all n since u(a,t) = (b,t) = 0
% lambda1 = alpha*dt/(2*dx), lambda2 = D*dt/dx^2.
%
% Reference book: D. Bradie, A friendly introduction to numerical analysis 
%
% Written by Tsogbayar Tsednee (PhD), Kyungpook National University, South Korea 
% Contact: tsog215@gmail.com
% Date: Jan 11, 2021
%
format short;
clear;
clc;
%%% ---
D = 1.00;
alpha = 1.0;
a = 0.; % x(0)
b = 1.0; % x(N+1)
N = 8;  % 
dx = (b-a)/N;
% grid and initial condition
% ---
x = zeros(N+1,1); % total number of points is N+1
for i = 1:N+1
    x(i) = a + (i-1)*dx;
end
%x;
% ---
ti = 0.; % t(0)
tf = 2.00; % t(N)
Nt = 32;
dt = (tf-ti)/Nt;
t = zeros(Nt,1); % total number of points is Nt
%
lambda1 = alpha*dt/(2.*dx);
lambda2 = D*dt/dx^2;
%
%  matrix equation is
% (I + lambda1*u_x + dt*diag(u) - lambda2*u_xx)*w^(n+1) = w^(n) + 
%                                                         lambda_1*b1^(n+1) + 
%                                                         lambda_2*b2^(n+1)
%                                                         dt*s^(n+1) 
%
%
u_mat_dx1 = zeros(N+1,N+1);
u_mat_dx2 = zeros(N+1,N+1);
for i = 2:N
    u_mat_dx2(i,i-1) = 1.0;
    u_mat_dx2(i,i) = -2.;
    u_mat_dx2(i,i+1) = 1.;
%
    u_mat_dx1(i,i-1) = -1.;
    u_mat_dx1(i,i) = 0.;
    u_mat_dx1(i,i+1) = 1.;
end
u_mat_dx2(1,1) = -2.;
u_mat_dx2(N+1,N+1) = -2. ;
u_mat_dx2(1,2) = 1.;
u_mat_dx2(N+1,N) = 1.;
%
u_mat_dx2 = u_mat_dx2(1:N+1,1:N+1);
u_mat_dx2 = lambda2*u_mat_dx2 ;
%
u_mat_dx1(1,1) = 0.;
u_mat_dx1(N+1,N+1) = 0.;
u_mat_dx1(1,2) = 1.;
u_mat_dx1(N+1,N) = -1.;
%
u_mat_dx1 = u_mat_dx1(1:N+1,1:N+1) ;
u_mat_dx1 = lambda1*u_mat_dx1;
%
unit_I = eye(N+1);
% at t0 = 0
w0 = x.^2; % % at t0 = 0.00 & initial
B_mat = 2.*unit_I;
evolution_mat = unit_I + u_mat_dx1 + dt.*B_mat - u_mat_dx2;
%%%
w0 = w0(2:N); % taking intou account the BC 

for j = 1:Nt
    t(j) = ti + (j-0)*dt;
    b_1n = [0.; zeros(N-3,1); -exp(-t(j))];    
    b_2n = [0.; zeros(N-3,1); exp(-t(j))];        
    s_n = exp(-t(j)).*(x(2:N).^2 + 2.*x(2:N) - 2.);   
%    
    w_old = w0 + lambda1.*b_1n + lambda2.*b_2n + dt.*s_n;
    w_new = evolution_mat(2:N,2:N)\(w_old  ); % % BTCS scheme
    w0 = w_new  ;
    %        
%    figure(1)
%    hold on
%    plot(x(2:N), w0, 'r-'); drawnow;    
%    plot(x(2:N), x(2:N).^2.*exp(-t(j)), 'bo'); drawnow;
%    hold off
    
%    
end
%%%
time = tf; % 
w_sol = [0.0; w_new; exp(-time)];
u_exact = x.^2.*exp(-time);
%
figure(1)
hold on
plot(x, w_sol, 'b')
plot(x, u_exact, 'ro')
hold off
%axis([0. 1. .00 1])
set(gca,'FontSize',18)
xlabel('x') % ,'fontsize',16
ylabel('u(x,2)') % ,'Rotation', 
%%%
%rms(w_sol - u_exact);


return
end
