function simulateLQG()
M = 1000;
m1 = 100;
m2 = 100;
l1 = 20;
l2 = 10;
g = 9.80;
A= [0 1 0 0 0 0 ; 0 0 -(m1*g/M) 0 -(m2*g/M) 0 ; 0 0 0 1 0 0 ; 0 0 -(g*(M+m1))/(M*l1) 0 -(m2*g)/(M*l1) 0 ;0 0 0 0 0 1; 0 0 -(m1*g)/(M*l2) 0 -(g*(M+m2))/(M*l2) 0]
B= [ 0; 1/M ;0; 1/(M*l1) ;0 ;1/(M*l2)];
C = [1 0 0 0 0 0];
D = 0;
Q=[1 0 0 0 0 0;0 1 0 0 0 0; 0 0 100 0 0 0; 0 0 0 1000 0 0; 0 0 0 0 150 0; 0 0 0 0 0 1500]
R=0.0001
K = lqr(A,B,Q,R);
sys_1 = ss(A,[B B],C,[zeros(1,1) zeros(1,1)]);
vd = 0.1 * eye(1);
vn = 1;
sen = [1];
known = [1];
[~,L,~] = kalman(sys_1,vd,vn,[],sen,known)
states = {'x','x_dot','theta1','theta1_dot','theta2','theta2_dot','e_1','e_2','e_3','e_4','e_5','e_6'};
inputs = {'F'};
outputs = {'x'};

Ac = [A-B*K B*K;zeros(size(A)) A-L*C];
Bc = zeros(12,1);
Cc = [C zeros(size(C))];
sys_cl_lqg = ss(Ac,Bc,Cc,D, 'statename',states,'inputname',inputs,'outputname',outputs);    
x0= [ 5 ; 0 ; 0.1 ; 0 ;0.2 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0]
t = 0:0.01:100;
F = zeros(size(t));
[Y,~,X] = lsim(sys_cl_lqg,F,t,x0);
figure
plot(t,Y(:,1),'b');

u = zeros(size(t));
for i = 1:size(X,1)
   u(i) = K * (X(i,1:6))';
end
Xhat = X(:,1) - X(:,6);
figure
subplot(3,1,1), plot(t,Xhat), hold on, plot(t,X(:,1),'r')

end
