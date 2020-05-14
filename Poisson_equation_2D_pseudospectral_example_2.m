function [] = Poisson_equation_2D_pseudospectral_example_2
%%% --- 2D Poisson problem with Dirichlet boundary condition
%
% u_xx + u_yy = (x^2+y^2)*cos(x*y) - cos(pi*x) for 0<x,y<1, y < 1, &
%                         u(x,0) = (cos(pi*x) - pi^2)/pi^2
%                         u(x,1) = (cos(pi*x) - pi^2*cos(x))/pi^2
%                         u(0,y) = (1/pi^2) - 1
%                         u(1,y) = -(1/pi^2 + cos(y)) 
%
% analytic solution: u(x,y) = (1/pi^2)*cos(pi*x) - cos(x*y)
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
% ----------------------------------------------------------------
clear; clc; format short
N = 32.; a = 1.; b = 0.;  % Initial computaional parameters; you may change them
%
[x,xw,D]=legDC2(N,a,b); % Differentiation matrix, nodes and weights 
%
D2 = D*D; % second-order differentiation matrix
%
y = x; % y-axis
%
[xx,yy] = meshgrid(x,y) ; % 2D-mesh
xx = xx(:); yy = yy(:);
D2 =(2/(b-a))^2*D2;
I = eye(N+1);
L = kron(I,D2) + kron(D2,I);

% Impose boundary conditions by replacing appropraite rows:
bxl = find(abs(xx) == 0 );
bxr = find(abs(xx) == 1 );
bub = find(abs(yy) == 0 );
but = find(abs(yy) == 1 );

%bc = find(abs(xx) == 0 | abs(yy) == 0| abs(xx) == 1 | abs(yy) == 1)      % boundary 
%b = find(abs(yy) ==1 | )
L(bxl,:) = zeros(size(bxl,1),(N+1)^2);
L(bxl,bxl) = eye(size(bxl,1));
%L;
%
L(bxr,:) = zeros(size(bxr,1),(N+1)^2);
L(bxr,bxr) = eye(size(bxr,1));
%L;
%
L(bub,:) = zeros(size(bub,1),(N+1)^2);
L(bub,bub) = eye(size(bub,1));
%L;
%
L(but,:) = zeros(size(but,1),(N+1)^2);
L(but,but) = eye(size(but,1));
%L;
%%% rhs
rhs = zeros((N+1)^2,1);
%
rhs(bxl) = (xx(bxl) == 0).*(cos(pi.*xx(bub)) - pi^2)/pi^2;
rhs(bxr) = (xx(bxr) == 1).*(cos(pi.*xx(but)) - pi^2*cos(xx(but)))/pi^2;
rhs(bub) = (yy(bub) == 0).*1./pi^2 -1.;
rhs(but) = (yy(but) == 1).*(-((1./pi^2) + cos(yy(bxr))));

%rhs = rhs

% Solve Laplace equation, reshape to @D, and plot:
%figure(1), clf, spy(L), drawnow
%L;
%u = L\rhs;
reshape(rhs,N+1,N+1) ;
%%%
[xxi,yyi] = meshgrid(x(2:N),y(2:N));
gs = zeros(N+1,N+1);
gs(2:N,2:N) = (xxi.^2 + yyi.^2).*cos(xxi.*yyi) - cos(pi.*xxi) ;
gs = reshape(gs,(N+1)*(N+1),1);
%
f = -gs + rhs;
u = L\f;
%rhs = reshape(rhs,N+1,N+1)
% Reshape long 1D resuls onto 2D grid:
uu = reshape(u,N+1,N+1)';
%%%
figure(1)
mesh(x,y,uu)  % plot the numerical solution
%
%%%
uexact = reshape((cos(pi*xx)./pi^2) - cos(xx.*yy),N+1,N+1);
%
figure(2)
rel_err = (uexact' - uu')./uexact'; % plot the relative error for solution
mesh(x,y,rel_err)

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



