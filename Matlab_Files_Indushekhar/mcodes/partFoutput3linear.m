%second Component (F) , Luenberger observer with output (x(t),theta_1(t),theta_2(t))%

m1=100;
m2=100;
M = 1000;
g = 9.8;
l1 =20;
l2 =10; 
A= [0 1 0 0 0 0 ; 0 0 -(m1*g/M) 0 -(m2*g/M) 0 ; 0 0 0 1 0 0 ; 0 0 -(g*(M+m1))/(M*l1) 0 -(m2*g)/(M*l1) 0 ;0 0 0 0 0 1; 0 0 -(m1*g)/(M*l2) 0 -(g*(M+m2))/(M*l2) 0]
B= [ 0; 1/M ;0; 1/(M*l1) ;0 ;1/(M*l2)];
C=[1 0 0 0 0 0 ; 0 0 1 0 0 0; 0 0 0 0 1 0]
D=[0]
x_0= [ 5 ; 0 ; 0.1 ; 0 ;0.2 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0]
 Q=[1 0 0 0 0 0;0 1 0 0 0 0; 0 0 100 0 0 0; 0 0 0 1000 0 0; 0 0 0 0 150 0; 0 0 0 0 0 1500]
R=0.0001
K=lqr(A,B,Q,R)
eig(A-B*K)
p=[-2 -3 -4 -5 -6 -7 ]
L=place(A',C',p)
L=L'
eig(A-L*C)
Ac=[(A-B*K) (B*K); zeros(size(A)) (A-L*C)];
Bc=[ B ;zeros(size(B))];
Cc= [C zeros(size(C))];
Dc=[0]
sys=ss(Ac,Bc,Cc,Dc)
initial(sys,x_0)
grid





