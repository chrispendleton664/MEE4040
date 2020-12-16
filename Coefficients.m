function Coefficients = GetCoefficients(Q, Qdot)

U0 = Q(1);
U1 = Qdot(1);


V0 = Q(2);
V1 = Qdot(2);


W0 = Q(3);
W1 = Qdot(3);


% Need to define qx, Pt, dt, Wn, U1, Un, V1, W1, qz

dt = 0.1;
a = 2;
b = 2;

% Blast Pressure Equation
pm = 3447E3;
t =  1;
tp =  0.1;
alpha = 2;
Pt = pm*(1-t/tp)*exp(-alpha*t/tp);
mdot = 0.1;
qz = 1;


% coefficients in the time dependent nonlinear differential equations note:
% d1, d2 and d3 are viscous damping coeffiients not sure if needed in this
% case or how to define 

a0 = a*b^9*mdot / 1260;

a1 = 3*a^2*Aijsum(3,3)*b^7+Aijsum(1,1)*b^9*pi^2/(315*a);

a2 = 9*a^4*b^4*(Aijsum(1,2)+Aijsum(3,3))/pi^6;

a3 = 16*b^3*(b^2*Bijsum(1,1)+a^2*(Bijsum(1,2)+2*Bijsum(3,3)))*(-12+pi^2)/(3*a^2*pi^3);

a4 = b^3*(Aijsum(1,1)*b^2*(45+pi^4)+a^2*(90*Aijsum(3,3)-Aijsum(1,2)*(-45+pi^4)))/(240*a^2*pi);

a5 = 0;

%b coefficients in the time dependent nonlinear differential equations 

b0 = 9*a^9*b*mdot/1260;

b1 = 3*a^7*b^2*Aijsum(3,3)+a^9*Aijsum(2,2)*pi^2/(315*b);

b2 = 9*a^4*b^4*(Aijsum(1,2)+Aijsum(3,3))/pi^6;

b3 = 16*a^3*(a^2*Bijsum(2,2)+b^2*(Bijsum(1,2)+2*Bijsum(3,3)))*(-12+pi^2)/(3*b^2*pi^3);

b4 = a^3*(Aijsum(2,2)*a^2*(45+pi^4)+b^2*(90*Aijsum(3,3)-Aijsum(1,2)*(-45+pi^4)))/(240*b^2*pi);

b5 = 0;

%c coefficients in the time dependent nonlinear differential equations 

c0 = a*b*mdot/4;

c1 = pi^4*(b^4*Dijsum(1,1)+a^4*Dijsum(2,2)+2*a^2*b^2*(Dijsum(1,2)+2*Dijsum(3,3)))/(4*a^3*b^3);

c2 = 8*pi^2*(Bijsum(1,2)-Bijsum(3,3))/(3*a*b);

c3 = (9*a^4*Aijsum(2,2)+2*a^2*(3*Aijsum(1,2)+4*Aijsum(3,3))*b^2+9*Aijsum(1,1)*b^4)*pi^4/(128*a^3*b^3);    

c4 = (16*b^3*(4*b^2*Bijsum(1,1)+a^2*(Bijsum(1,2)+2*Bijsum(3,3))))*(-12+pi^2)/(3*a^2*pi^3);

c5 = (16*a^3*(4*a^2*Bijsum(3,3)+b^2*(Bijsum(1,2)+2*Bijsum(3,3))))*(-12+pi^2)/(3*b^2*pi^3);

c6 = (3*b^3*(Aijsum(1,2)+Aijsum(3,3)))/(4*pi)  +  (Aijsum(1,1)*((b^5/60)+(3*b^5/(4*pi^4)))*pi^3)/(2*a^2)  -  (Aijsum(1,2)*(b^5/60+3*b^5/(4*pi^4))*pi^3)/(2*b^2);  

c7 = (3*a^3*(Aijsum(1,2)+Aijsum(3,3)))/(4*pi)  +  (Aijsum(2,2)*((5^5/60)+(3*a^5/(4*pi^4)))*pi^3)/(2*b^2)  -  (Aijsum(1,2)*(a^5/60+3*a^5/(4*pi^4))*pi^3)/(2*a^2);

c8 = -4*a*b*qz/pi^2;

% Coefficients in the finite difference equations

Coefficients.A1 = a0/dt^2+a1;
Coefficients.A2 = a2;
Coefficients.A3 = a3+a4*W0;
Coefficients.A4 = a0*U1/dt+a0*U0/dt^2-a5;

Coefficients.B1 = b2;
Coefficients.B2 = b0/dt^2+b1;
Coefficients.B3 = b3+b4*W0;
Coefficients.B4 = b0*V1/dt+b0*V0/dt^2-b5;

Coefficients.C1 = c4+c6*W0;
Coefficients.C2 = c5+c7*W0;
Coefficients.C3 = c0/dt^2+c1+c2*W0+c3*W0^2;
Coefficients.C4 = c0*W1/dt+c0*W0/dt^2+c8;

Coefficients.D1 = B2-B1*A2/A1;
Coefficients.D2 = B3-B1*A3/A1;
Coefficients.D3 = B4-B1*A4/A1;

Coefficients.E11 = C2-C1*A2/A1;
Coefficients.E22 = C3-C1*A3/A1;
Coefficients.E3 = C4-C1*A4/A1;
Coefficients.E4 = E11*D3/D1;
Coefficients.E5 = E22- E11*D2/D1;
end