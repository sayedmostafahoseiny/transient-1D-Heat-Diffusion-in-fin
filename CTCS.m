clear all ; close all ; clc ;

% definition of parameters %
L = 0.2 ; % length of fin
R = 0.01 ; % radius of fin
h = 25 ; % convection heat transfer coefficient of ambient
C = 897 ; % specific heat capacity of fin
rho = 2700 ; % density of fin
k = 200 ; % conduction heat transfer coefficient
Tinf = 300 ; % ambient temperature in kelvin
Tbase = 500 ; % temperature of base of fin
teta_base = 1 ; % dimensionless temperature at base
P = 2*pi*R ; % perimeter of fin
A_s = pi*R^2 ; % cross section area of fin
m = (h*P/(k*A_s))^0.5 ;
% mesh creation
n = 10 ; % number of divisions
deltax = 1/n ;
alpha = k/(rho*C);
nu = 0.5 ;
deltat = nu*(deltax^2) ;
NTS = 1000; % number of time steps

B = deltat/(deltax^2);
S =(1-2*B-deltat*(m*L)^2);

% initial condition
theta(1,1)=1; % for left boundary node
for j=2:n+1
theta(1,j)=0 ;
end

t = 0 ;
% for using CTCS we have to march first step by BTCS or FTCS, we use BTCS here
% assembling coefficient matrix A in sparse form
A=sparse([],[],[],n+1,n+1); % allocation of sparse matrix A
% time loop
t = t + deltat ;
A(1,1)=1;% left boundary node
b(1,1)= 1 ;% left boundary condition
for i=2:n ; % interior nodes
    A(i,i)= S -2 ;
    A(i,i-1)= B;
    A(i,i+1)= B ;
    b(i,1)= -theta(1,i) ;
end
b(n+1,1) = -theta(1,n+1); % right boundary condition
% right boundary node
A(n+1,n)= 2*B ;
A(n+1,n+1)= S-2-(2*B*h*L*deltax)/k ;

% solving system of equations
% solution=gmres(A,b,[],1e-7,1000);
solution=A\b;
theta(2,:)= solution' ;
theta=theta';

% we march other steps by CTCS discretization method
U=sparse([],[],[],n+1,n+1); % allocation of sparse matrix A
V=sparse([],[],[],n+1,n+1);
U(1,1)=1;% left boundary node
for i=2:n ; % interior nodes
    U(i,i)= -4*B-2*deltat*((m*L)^2);
    U(i,i-1)= 2*B;
    U(i,i+1)= 2*B ;
    V(i,i) = 1 ;
end
% right boundary node
U(n+1,n)= 4*B ;
U(n+1,n+1)= -((4*B*h*L*deltax/k) + 4*B + 2*deltat*((m*L)^2)) ;
V(n+1,n+1) = 1 ;
% time loop
for j=2:NTS % "j" stands for time and "i" stands for space
t = t + deltat ;
theta(:,j+1)= U * theta(:,j) + V * theta (:,j-1) ; % left boundary node
 if sum(abs(theta(:,j+1)-theta(:,j))) < 0.00001 % criterion for reaching to steady state
    disp ( ' solution has reached to steady state after ' )
    t
    %break
end
end

% postprocessing
x=linspace(0,1,n+1);
plot(x,theta(:,end),'b*','linewidth',0.5)
grid
xlabel(' x* ' )
ylabel(' theta ' )
