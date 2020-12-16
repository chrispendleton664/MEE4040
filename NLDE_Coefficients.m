% Refer to main file

addpath(AttemptAtInstalling3DArrays)

% Need to define qx, Pt, dt, Wn, Udotn, Un, Vdotn, Wdotn,

% Blast Pressure Equation
pm =  %Peak Blast Pressure
t =  %time
tp =  %positive phase duration
alpha = %waveform parameter
Pt = pm*(1-t/tp)*e^(-alpha*t/tp);


% coefficients in the time dependent nonlinear differential equations note:
% d1, d2 and d3 are viscous damping coeffiients not sure if needed in this
% case or how to define 

a0 = a^11*b*mdot / 18480;

adash = a^11*b*d1/18480;

a1 = a^9*(33*Aijsum(1,1)*b^2 + a^2*Aijsum(6,6)*pi^2)/(13680*b);

a2 = a^5*(Aijsum(1,2)+Aijsum(6,6))*b^5*(-15 + pi^2)^2/(4*pi^8);

a3 = -a^3*(3*b^2*Bijsum(1,1) + a^2*(Bijsum(1,2)+2*Bijsum(6,6)))*(-15+pi^2)/b*pi^2;

a4 = 5*a^3*(b^2*Aijsum(2,2)*(15-4*pi^2)+3*(Aijsum(1,2)-Aijsum(6,6))*a^2*(63-4*pi^2))/64*b*pi^2;

a5 = -a*b*qx;

%b coefficients in the time dependent nonlinear differential equations 

b0 = a*b^11*mdot/18480;

bdash = a^11*b*d2/18480;

b1 = b^9*(33*Aijsum(2,2)*a^2 + b^2*Aijsum(6,6)*pi^2)/(13680*b);

b2 = a^5*(Aijsum(1,2)+Aijsum(6,6))*b^5*(-15 + pi^2)^2/(4*pi^8);

b3 = b^3*(3*a^2*Bijsum(2,2) + b^2*(Bijsum(1,2)+2*Bijsum(6,6)))*(-15+pi^2)/a*pi^2;

b4 = 5*b^3*(a^2*Aijsum(2,2)*(15-4*pi^2)+3*(Aijsum(1,2)-Aijsum(6,6))*b^2*(63-4*pi^2))/64*a*pi^2;

b5 = -a*b*qy;

%c coefficients in the time dependent nonlinear differential equations 

c0 = 9*a*b*mdot/4;

cdash = 9*a*b*d3/4;

c1 = 4*pi^4*(3*b^4*Dijsum(1,1)+3*a^4*Dijsum(2,2)+2*a^2*b^2*(Dijsum(1,2)+2*Dijsum(6,6)))/a^3*b^3;

c2 = 24*pi^4*(Bijsum(1,2)-Bijsum(6,6))/a*b;

c3 = 20*pi^4*(21*a^4*Aijsum(2,2)+10*a^2*(3*Aijsum(1,2)+4*Aijsum(6,6))*b^2+21*Aijsum(1,1)*b^4)/32*a^3*b^3;    

c4 = -a^3*(3*b^2*Bijsum(1,1)+a^2*(Bijsum(1,2)+2*Bijsum(6,6)))*(-15+pi^2)/b*pi^2;

c5 = -b^3*(3*a^2*Bijsum(2,2)+b^2*(Bijsum(1,2)+2*Bijsum(6,6)))*(-15+pi^2)/a*pi^2;

c6 = 5*a^3*(Aijsum(1,1)*b^2*(15-4*pi^2)-3*a^2*(Aijsum(1,2)-Aijsum(6,6))*(-63+4*pi^2))/32*b*pi^2;

c7 = 5*b^3*(Aijsum(2,2)*a^2*(15-4*pi^2)-3*b^2*(Aijsum(1,2)-Aijsum(6,6))*(-63+4*pi^2))/32*a*pi^2;

c8 = -a*b*Pt;

% Coefficients in the finite difference equations

A1 = a0/dt^2+a1+adash/dt;
A2 = a2;
A3 = a3+a4*Wn;
A4 = a0*Udotn/dt+a0*Un/dt^2+adash*Un/dt-a5;

B1 = b2;
B2 = b0/dt^2+bdash/dt+b1;
B3 = b3+b4*Wn;
B4 = b0*Vdotn/dt+b0*Vn/dt^2+bdash*Vn/dt-b5;

C1 = c4+c6*Wn;
C2 = c5+c7*Wn;
C3 = c0/dt^2+c1+c2*Wn+c3*Wn^2+cdash/dt;
C4 = c0*Wdotn/dt+(c0/dt^2+cdash/dt)*Wn+c8;

D1 = B2-B1*A2/A1;
D2 = B3-B1*A3/A1;
D3 = B4-B1*A4/A1;

E1 = C2-C1*A2/A1;
E2 = C3-C1*A3/A1;
E3 = C4-C1*A4/A1;
E4 = E1*D3/D1;
E5 = E2- E1*D2/D1;


% Different shape (approximation) functions 

% Uo = 

% Stating the Linear and Nonlinear Coefficients for the Equations of Motion

% L11 = -Aijsum(1,1)