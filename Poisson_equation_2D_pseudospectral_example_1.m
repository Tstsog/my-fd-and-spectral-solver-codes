function [] = Poisson_equation_2D_pseudospectral_example_1
%%% --- 2D Poisson problem with Dirichlet boundary condition
%
% u_xx + u_yy = 10*sin(8*x*(y-1)) for -1 < x, y < 1, u = 0 on the baundary
% analytic solution: 
%
% The numerical method: pseudospectral method based on Legendre-spectral
% method.
% Uses: function legDC2 - differentation matrix based on
% Legendre-Gauss-Lobatto nodes
% Reference book: N.Trefethen, Spectral Methods in MATLAB SIAM, Jan 1, 2000 & 
% Tsogbayar Tsednee, PhD thesis, available at https://yorkspace.library.yorku.ca/xmlui/handle/10315/28160  
%
% Written by Tsogbayar Tsednee (PhD), California State University Northridge 
% Contact: tsog215@gmail.com
% Date: August 15, 2012 & May 12, 2020
%%%
% ----------------------------------------------------------------
clear; clc; format long
N = 64.; a = -1.; b = 1.;  % Initial computaional parameters; you may change them
%
[x,xw,D]=legDC2(N,a,b); % Differentiation matrix, nodes and weights 
%
D2 = D*D; % second-order differentiation matrix
%
y = x; % y-axis
%
[xx,yy] = meshgrid(x(2:N),y(2:N)); % 2D-mesh
xx = xx(:); yy = yy(:);
% right-hand side part
f = 10*sin(8*xx.*(yy-1)); 
%
D2 =(2/(b-a))^2*D2;
D2 = D2(2:N,2:N); % taking into account the boundary condition
I = eye(N-1);     % unit matrix
L = kron(I,D2) + kron(D2,I);      % 2D-Laplacian with differentation matrix
%
%figure(1), clf, spy(L), drawnow
tic, u = L\f; toc                 % solve problem and watch the clock
%u;
% Reshape long 1D resuls onto 2D grid:
uu = zeros(N+1,N+1);
uu(2:N,2:N) = reshape(u,N-1,N-1); 
%
figure(2)
mesh(x,y,uu)  % plot the numerical solution

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xi,w,D]=legDC2(N,a,b)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% legDc.m
%
% Computes the Legendre differentiation matrix with collocation at the 
% Legendre-Gauss-Lobatto nodes.
%
% Reference: 
%   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
%   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
%
% Written by Greg von Winckel - 05/26/2004
% Contact: gregvw@chtm.unm.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Truncation + 1
N1=N+1;

% CGL nodes
xc=cos(pi*(0:N)/N)';

% Uniform nodes
xu=linspace(-1,1,N1)';

% Make a close first guess to reduce iterations
if N<3
    x=xc;
else
    x=xc+sin(pi*xu)./(4*N);
end

% The Legendre Vandermonde Matrix
P=zeros(N1,N1);

% Compute P_(N) using the recursion relation
% Compute its first and second derivatives and 
% update x using the Newton-Raphson method.

xold=2;
while max(abs(x-xold))>eps

    xold=x;
        
    P(:,1)=1;    P(:,2)=x;
    
    for k=2:N
        P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
    end
     
    x=xold-( x.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );
end

X=repmat(x,1,N1);
Xdiff=X-X'+eye(N1);

L=repmat(P(:,N1),1,N1);
L(1:(N1+1):N1*N1)=1;
D=(L./(Xdiff.*L'));
D(1:(N1+1):N1*N1)=0;
D(1)=(N1*N)/4;
D(N1*N1)=-(N1*N)/4;

% Linear map from[-1,1] to [a,b]
xi=(a*(1-x)+b*(1+x))/2;      % added by Tsogbayar Tsednee

% Compute the weights
w=(b-a)./(N*N1*P(:,N1).^2); % added by Tsogbayar Tsednee
return
end

