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
n = 640 ; % number of divisions
deltax = 1/n ;
alpha = k/(rho*C);
nu = 0.5 ;
deltat = 0.001 ;
NTS = 10000; % number of time steps

B = deltat/(deltax^2);
S =(1-2*B-deltat*(m*L)^2);
E =2*B;
F = S-(2*h*L*deltat)/(k*deltax);

% initial condition
%theta(1,1)=1; % for left boundary node
for j=1:n+1
theta(1,j)=0 ;
end

% assembling coefficient matrix A in sparse form
A=sparse([],[],[],n+1,n+1); % allocation of sparse matrix A
% time loop
t = 0 ;
for j=1:NTS % "j" stands for time and "i" stands for space
t = t + deltat ;

A(1,1)=1;% left boundary node
b(1,1)= 1 ;% left boundary condition
for i=2:n ; % interior nodes
    A(i,i)= S -2 ;
    A(i,i-1)= B;
    A(i,i+1)= B ;
    b(i,1)= -theta(j,i) ;
end
b(n+1,1) = -theta(j,n+1); % right boundary condition
% right boundary node
A(n+1,n)= 2*B ;
A(n+1,n+1)= S-2-(2*B*h*L*deltax)/k ;

% solving system of equations
% solution=gmres(A,b,[],1e-7,1000);
solution=A\b;
theta(j+1,:)= solution' ;

if sum(abs(theta(j+1,:)-theta(j,:))) < 0.00001 % condition for reaching to steady state
    disp ( ' solution has reached to steady state after ' )
    t
    break
end
end

% postprocessing
x=linspace(0,1,n+1);
plot(x,theta(end,:),'b*','linewidth',0.5);
grid
xlabel(' x* ' )
ylabel(' theta ' )

% analytical solution
xx=linspace(0,L,length(x));
teta_Analytic=(cosh(m*(L-xx))+(h/(m*k))*sinh(m*(L-xx)))/(cosh(m*L)+(h/(m*k))*sinh(m*L));
figure(2)
plot1=plot(xx,theta(end,:),'b*','linewidth',0.5); % steady state solution
hold on
plot2=plot(xx,teta_Analytic,'r','linewidth',2.5);
grid on
title( ' Temperature distribution ' )
xlabel( ' x (m) ' )
ylabel( ' theta ' )
legend ([plot1 plot2],{'numerical (BTCS)','Analitical solution'}) % validation

% Error Analysis
L1_norm = sum(abs(teta_Analytic-theta(end,:)))/(n+1) % calculation of first norm of Error
L2_norm = sqrt(sum((teta_Analytic-theta(end,:)).^2)/(n+1)) % calculation of second norm of Error

% post processing
T=theta(end,:)*(Tbase-Tinf)+Tinf; % temperature distribution in (C) in fin
q=-k*((T(2)-T(1))/(xx(2)-xx(1))) % the flux which enter the fin from base
eta=(pi*q*R^2)/(h*2*pi*R*L*(Tbase-Tinf))*100 % fin efficiency
eps = (pi*q*R^2)/(pi*R^2*h*(Tbase-Tinf)) % fin performance coefficient
T_ave=sum(T)/length(T)
