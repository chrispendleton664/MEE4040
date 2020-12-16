% Refer to main file

addpath(AttemptAtInstalling3DArrays)

function F = root2d(t)

syms U V W;
U0 = U;
U1 = diff(U0, 1);
U2 = diff(U0, 2);

V0 = V;
V1 = diff(V0, 1);
V2 = diff(V0, 2);

W0 = W;
W1 = diff(W0, 1);
W2 = diff(W0, 2);


a= 2.54;
b= 2.54;
h=0.17;



% Need to define qx, Pt, dt, Wn, U1, Un, V1, W1, qz

% Blast Pressure Equation
pm = 3447E3;
t =  1;
tp =  0.1;
alpha = 2;
Pt = pm*(1-t/tp)*e^(-alpha*t/tp)
mdot = 1;
qz = 1;


% coefficients in the time dependent nonlinear differential equations note:
% d1, d2 and d3 are viscous damping coeffiients not sure if needed in this
% case or how to define 

a0 = a*b^9*mdot / 1260;

a1 = 3*a^2*Aijsum(6,6)*b^7+Aijsum(1,1)*b^9*pi^2/(315*a);

a2 = 9*a^4*b^4*(Aijsum(1,2)+Aijsum(6,6))/pi^6;

a3 = 16*b^3*(b^2*Bijsum(1,1)+a^2*(Bijsum(1,2)+2*Bijsum(6,6)))*(-12+pi^2)/(3*a^2*pi^3);

a4 = b^3*(Aijsum(1,1)*b^2*(45+pi^4)+a^2*(90*Aijsum(6,6)-Aijsum(1,2)*(-45+pi^4)))/(240*a^2*pi);

a5 = 0;

%b coefficients in the time dependent nonlinear differential equations 

b0 = 9*a^9*b*mdot/1260;

b1 = 3*a^7*b^2*Aijsum(6,6)+a^9*Aijsum(2,2)*pi^2/(315*b);

b2 = 9*a^4*b^4*(Aijsum(1,2)+Aijsum(6,6))/pi^6;

b3 = 16*a^3*(a^2*Bijsum(2,2)+b^2*(Bijsum(1,2)+2*Bijsum(6,6)))*(-12+pi^2)/(3*b^2*pi^3);

b4 = a^3*(Aijsum(2,2)*a^2*(45+pi^4)+b^2*(90*Aijsum(6,6)-Aijsum(1,2)*(-45+pi^4)))/(240*b^2*pi);

b5 = 0;

%c coefficients in the time dependent nonlinear differential equations 

c0 = a*b*mdot/4;

c1 = pi^4*(b^4*Dijsum(1,1)+A^4*Dijsum(2,2)+2*a^2*b^2*(Dijsum(1,2)+2*Dijsum(6,6)))/(4*a^3*b^3);

c2 = 8*pi^2*(Bijsum(1,2)-Bijsum(6,6))/(3*a*b);

c3 = (9*a^4*Aijsum(2,2)+2*a^2*(3*Aijsum(1,2)+4*Aijsum(6,6))*b^2+9*Aijsum(1,1)*b^4)*pi^4/(128*a^3*b^3);    

c4 = (16*b^3*(4*b^2*Bijsum(1,1)+a^2*(Bijsum(1,2)+2*Bijsum(6,6))))*(-12+pi^2)/(3*a^2*pi^3);

c5 = (16*a^3*(4*a^2*Bijsum(2,2)+b^2*(Bijsum(1,2)+2*Bijsum(6,6))))*(-12+pi^2)/(3*b^2*pi^3);

c6 = (3*b^3*(Aijsum(1,2)+Aijsum(6,6)))/(4*pi)  +  (Aijsum(1,1)*((b^5/60)+(3*b^5/(4*pi^4)))*pi^3)/(2*a^2)  -  (Aijsum(1,2)*(b^5/60+3*b^5/(4*pi^4))*pi^3)/(2*b^2);  

c7 = (3*a^3*(Aijsum(1,2)+Aijsum(6,6)))/(4*pi)  +  (Aijsum(2,2)*((5^5/60)+(3*a^5/(4*pi^4)))*pi^3)/(2*b^2)  -  (Aijsum(1,2)*(a^5/60+3*a^5/(4*pi^4))*pi^3)/(2*a^2);

c8 = -4*a*b*qz/pi^2;

% Coefficients in the finite difference equations

A1 = a0/dt^2+a1;
A2 = a2;
A3 = a3+a4*Wn;
A4 = a0*U1/dt+a0*Un/dt^2-a5;

B1 = b2;
B2 = b0/dt^2+b1;
B3 = b3+b4*Wn;
B4 = b0*V1/dt+b0*Vn/dt^2-b5;

C1 = c4+c6*Wn;
C2 = c5+c7*Wn;
C3 = c0/dt^2+c1+c2*Wn+c3*Wn^2;
C4 = c0*W1/dt+c0*Wn/dt^2+c8;

D1 = B2-B1*A2/A1;
D2 = B3-B1*A3/A1;
D3 = B4-B1*A4/A1;

E1 = C2-C1*A2/A1;
E2 = C3-C1*A3/A1;
E3 = C4-C1*A4/A1;
E4 = E1*D3/D1;
E5 = E2- E1*D2/D1;


F(1) = a0*U2 + a1*U0 + a2*V0 + a3*W0 + a4*W0^2 + a5 = 0;
F(2) = b0*V2 + b1*V0 + b2*U0 + b3*W0 + b4*W0^2 + b5 = 0;
F(3) = c0*W2 + c1*W0 + c2W0^2 + c3*W0^3 + c4*U0 + c5*V0 + c6*U0*W0 + c8 = 0;
fun = @root2d;
k0 = [0,0];
k = fsolve(fun,k0)
end




% Different shape (approximation) functions 

% Uo = 

% Stating the Linear and Nonlinear Coefficients for the Equations of Motion

% L11 = -Aijsum(1,1)