clear all
close all
clc
tic

% =========================================================
syms th1(t) th2(t) l3(t) th4(t) th5(t) th6(t) t
syms th1p(t) th2p(t) l3p(t) th4p(t) th5p(t) th6p(t)

syms m1 m2 m3 m4 m5 m6
syms Ixx1 Iyy1 Izz1 Ixx2 Iyy2 Izz2 Ixx3 Iyy3 Izz3
syms Ixx4 Iyy4 Izz4 Ixx5 Iyy5 Izz5 Ixx6 Iyy6 Izz6

syms l1 l2 l4 l5 l6 lc1 lc2 lc3 lc4 lc5 lc6
syms g pi

Qp = [th1p; th2p; l3p; th4p; th5p; th6p];
Qp = Qp(t);


RP = [0 0 1 0 0 0];

% Cinemática directa

P(:,:,1)=[0;0;l1];
R(:,:,1)=[cos(th1) -sin(th1) 0;
          sin(th1)  cos(th1) 0;
          0 0 1]*[1 0 0;0 0 -1;0 1 0];

P(:,:,2)=[0;0;-l2];
R(:,:,2)=[cos(th2) -sin(th2) 0;
          sin(th2)  cos(th2) 0;
          0 0 1]*[1 0 0;0 0 1;0 -1 0];

P(:,:,3)=[0;0;l3];
R(:,:,3)=eye(3);

P(:,:,4)=[0;0;l4];
R(:,:,4)=[cos(th4) -sin(th4) 0;
          sin(th4)  cos(th4) 0;
          0 0 1]*[1 0 0;0 0 -1;0 1 0];

P(:,:,5)=[-l5*cos(th5); l5*sin(th5); 0];
R(:,:,5)=[cos(th5) -sin(th5) 0;
          sin(th5)  cos(th5) 0;
          0 0 1]*[1 0 0;0 0 1;0 -1 0];

P(:,:,6)=[0;0;l6];
R(:,:,6)=[cos(th6) -sin(th6) 0;
          sin(th6)  cos(th6) 0;
          0 0 1];

Z=zeros(1,3);

T(:,:,1)=[R(:,:,1) P(:,:,1);Z 1];
for i=2:6
    T(:,:,i)=T(:,:,i-1)*[R(:,:,i) P(:,:,i);Z 1];
end

for i=1:6
    RO(:,:,i)=T(1:3,1:3,i);
    PO(:,:,i)=T(1:3,4,i);
end

% ======== ESLABÓN 1 ========
Jv1=cross([0;0;1],PO(:,:,1));
Jw1=[0;0;1];
V1=Jv1*Qp(1);
W1=Jw1*Qp(1);

disp('======  1 ======')
disp('Velocidad lineal:')
pretty(V1)
disp('Velocidad angular:')
pretty(W1)

P01=subs(P(:,:,1)/2,l1,lc1);
I1=diag([Ixx1 Iyy1 Izz1]);
K1=1/2*m1*(V1+cross(W1,P01)).'*(V1+cross(W1,P01)) + 1/2*W1.'*I1*W1;
disp('Energía cinética:')
pretty(K1)

% ======== ESLABÓN 2 ========
Jv2=[cross([0;0;1],PO(:,:,2)) cross(RO(:,3,1),PO(:,:,2)-PO(:,:,1))];
Jw2=[[0;0;1] RO(:,3,1)];
V2=Jv2*Qp(1:2);
W2=Jw2*Qp(1:2);

disp('======  2 ======')
disp('Velocidad lineal:')
pretty(V2)
disp('Velocidad angular:')
pretty(W2)

P12=subs(P(:,:,2)/2,l2,lc2);
I2=diag([Ixx2 Iyy2 Izz2]);
K2=1/2*m2*(V2+cross(W2,P12)).'*(V2+cross(W2,P12)) + 1/2*W2.'*I2*W2;
disp('Energía cinética:')
pretty(K2)

% ======== ESLABÓN 3 ========
Jv3=[Jv2 RO(:,3,2)];
Jw3=[Jw2 [0;0;0]];
V3=Jv3*Qp(1:3);
W3=Jw3*Qp(1:3);

disp('======  3 ======')
disp('Velocidad lineal:')
pretty(V3)
disp('Velocidad angular:')
pretty(W3)

P23=subs(P(:,:,3)/2,l3,lc3);
I3=diag([Ixx3 Iyy3 Izz3]);
K3=1/2*m3*(V3+cross(W3,P23)).'*(V3+cross(W3,P23)) + 1/2*W3.'*I3*W3;
disp('Energía cinética:')
pretty(K3)

% ======== ESLABÓN 4 ========
Jv4=[Jv3 cross(RO(:,3,3),PO(:,:,4)-PO(:,:,3))];
Jw4=[Jw3 RO(:,3,3)];
V4=Jv4*Qp(1:4);
W4=Jw4*Qp(1:4);

disp('======  4 ======')
disp('Velocidad lineal:')
pretty(V4)
disp('Velocidad angular:')
pretty(W4)

P34=subs(P(:,:,4)/2,l4,lc4);
I4=diag([Ixx4 Iyy4 Izz4]);
K4=1/2*m4*(V4+cross(W4,P34)).'*(V4+cross(W4,P34)) + 1/2*W4.'*I4*W4;
disp('Energía cinética:')
pretty(K4)

% ======== ESLABÓN 5 ========
Jv5=[Jv4 cross(RO(:,3,4),PO(:,:,5)-PO(:,:,4))];
Jw5=[Jw4 RO(:,3,4)];
V5=Jv5*Qp(1:5);
W5=Jw5*Qp(1:5);

disp('======  5 ======')
disp('Velocidad lineal:')
pretty(V5)
disp('Velocidad angular:')
pretty(W5)

P45=subs(P(:,:,5)/2,l5,lc5);
I5=diag([Ixx5 Iyy5 Izz5]);
K5=1/2*m5*(V5+cross(W5,P45)).'*(V5+cross(W5,P45)) + 1/2*W5.'*I5*W5;
disp('Energía cinética:')
pretty(K5)

% ======== ESLABÓN 6 ========
Jv6=[Jv5 cross(RO(:,3,5),PO(:,:,6)-PO(:,:,5))];
Jw6=[Jw5 RO(:,3,5)];
V6=Jv6*Qp(1:6);
W6=Jw6*Qp(1:6);

disp('====== ESLABÓN 6 ======')
disp('Velocidad lineal:')
pretty(V6)
disp('Velocidad angular:')
pretty(W6)

P56=subs(P(:,:,6)/2,l6,lc6);
I6=diag([Ixx6 Iyy6 Izz6]);
K6=1/2*m6*(V6+cross(W6,P56)).'*(V6+cross(W6,P56)) + 1/2*W6.'*I6*W6;
disp('Energía cinética:')
pretty(K6)

K_Total = simplify(K1 + K2 + K3 + K4 + K5 + K6);
disp('Energía Cinética Total');
pretty(K_Total);

% Energia potencial p = m g h

h1 = P01(2);
h2 = P12(2);
h3 = P23(2);
h4 = P34(2);
h5 = P45(2);
h6 = P56(2);

U1 = m1*g*h1;
U2 = m2*g*h2;
U3 = m3*g*h3;
U4 = m4*g*h4;
U5 = m5*g*h5;
U6 = m6*g*h6;

% Energía potencial total
U_Total = U1 + U2 + U3 + U4 + U5 + U6;

disp('Energía Potencial Total');
pretty(U_Total);

% lagranjiano
Lagrangiano = simplify(K_Total - U_Total);
disp('Lagrangiano');
pretty(Lagrangiano);

% ================= ENERGÍA TOTAL =================
H = simplify(K_Total + U_Total);
disp('Energía Total');
pretty(H);