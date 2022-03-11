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
nu = 0.5 ; % ( CFL condition )
deltat = (nu*deltax^2); % calculating time step based on CFL condition
NTS = 1000000; % number of time steps

B = deltat/(deltax^2);
S =(1-2*B-deltat*((m*L)^2));
E =2*B;
F = S-((2*h*L*deltat)/(k*deltax));

% initial condition
theta(1,1)=1; % for left boundary node
for i=2:n+1
theta(i,1)=0 ;
end

A=sparse([],[],[],n+1,n+1); % allocation of sparse matrix A
A(1,1)=1;% left boundary node
for i=2:n ; % interior nodes
    A(i,i)= S ;
    A(i,i-1)= B;
    A(i,i+1)= B ;
end
% right boundary node
A(n+1,n)= E ;
A(n+1,n+1)= F ;
% time loop
t = 0 ;
for j=1:NTS % "j" stands for time and "i" stands for space
t = t + deltat ;
theta(:,j+1)= A * theta(:,j) ; % left boundary node
 if sum(abs(theta(:,j+1)-theta(:,j))) < 0.00001 % criterion for reaching to steady state
    disp ( ' solution has reached to steady state after ' )
    t
    break
end
end

x=linspace(0,1,n+1);
plot(x',theta(:,end),'b*','linewidth',1)
grid
xlabel( ' x* ')
ylabel( ' theta ' )


% analytical solution
xx=linspace(0,L,length(x));
teta_Analytic=(cosh(m*(L-xx))+(h/(m*k))*sinh(m*(L-xx)))/(cosh(m*L)+(h/(m*k))*sinh(m*L));
figure(2)
plot1=plot(xx,theta(:,end),'b*','linewidth',0.5); % steady state solution
hold on
plot2=plot(xx,teta_Analytic,'r','linewidth',2.5);
grid on
title( ' Temperature distribution ' )
xlabel( ' x (m) ' )
ylabel( ' theta ' )
legend ([plot1 plot2],{'numerical (FTCS)','Analitical solution'}) % validation
