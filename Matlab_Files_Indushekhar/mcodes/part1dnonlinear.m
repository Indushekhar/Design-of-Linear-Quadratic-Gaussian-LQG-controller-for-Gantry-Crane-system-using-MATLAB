%Component1 D , Non linear response to initial condition#

m1=100;
m2=100;
M = 1000;
g = 9.8;
l1 =20;
l2 =10; 
x_0= [ 5 ; 0; 0.1 ;  0 ; 0.2; 0]
t=0:0.01:100;%Timestep & Final Time
[t,x]=ode45(@findx,t,x_0);
plot(t,x,'linewidth',1.5);
title('Nonlinear Response');
function dx= findx(t,x)
m1=100;
m2=100;
M = 1000;
g = 9.8;
l1 =20;
l2 =10; 
A= [0 1 0 0 0 0 ; 0 0 -(m1*g/M) 0 -(m2*g/M) 0 ; 0 0 0 1 0 0 ; 0 0 -(g*(M+m1))/(M*l1) 0 -(m2*g)/(M*l1) 0 ;0 0 0 0 0 1; 0 0 -(m1*g)/(M*l2) 0 -(g*(M+m2))/(M*l2) 0];
B= [ 0; 1/M ;0; 1/(M*l1) ;0 ;1/(M*l2)];
C=[1 0 0 0 0 0];
D=[0];
Q=[1 0 0 0 0 0;0 1 0 0 0 0; 0 0 100 0 0 0; 0 0 0 1000 0 0; 0 0 0 0 150 0; 0 0 0 0 0 1500]
R=0.0001;
K=lqr(A,B,Q,R);
eig(A-B*K);
u = -K*x;
dx = zeros(6,1);
dx(1)= x(2);
dx(2)= -((-u) + m1*l1*x(4)^2*sin(x(3)) + m1*g*sin(x(3))*cos(x(3))+ m2*l2*x(6)^2*sin(x(5))+m2*g*sin(x(5))*cos(x(5)))/(M+m1*sin(x(3)^2)+ m2*sin(x(5)^2));
dx(3)= x(4);
dx(4)= -((-u) +(M+m1)*g*sin(x(3)) + m1*l1*x(4)^2*sin(x(3))*cos(x(3)) + m2*l2*x(6)^2*sin(x(5))*cos(x(3)) + m2*g*sin(x(5))*cos(x(3)-x(5)))/(( M+ m1*sin(x(3)^2)+ m2*sin(x(5)^2))*l1);
dx(5)= x(6);
dx(6)= -((-u) + m1*l1*x(4)^2*sin(x(3))*cos(x(5)) +m1*g*sin(x(3))*cos(x(3)-x(5)) + (M+m1)*g*sin(x(5))+m2*l2*x(6)^2*sin(x(5))*cos(x(5)))/(( M+ m1*sin(x(3)^2)+ m2*sin(x(5)^2))*l2);
end

