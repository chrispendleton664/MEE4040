fde = struct(A1, A2, A3, A4, B1, B2, B3, B4, C1, C2, C3, C4, D1, D2, D3, E1, E2, E3, E4, E5)

% Coefficients in the finite difference equations
sys (t) 

A1 = a0/dt^2+a1;
A2 = a2;
A3 = a3+a4*Wn;
A4 = a0*Un*d/dt*d/dt+a0*Un/dt^2-a5;

B1 = b2;
B2 = b0/dt^2+b1;
B3 = b3+b4*Wn;
B4 = b0*V*d/dt*d/dt+b0*Vn/dt^2-b5;

C1 = c4+c6*Wn;
C2 = c5+c7*Wn;
C3 = c0/dt^2+c1+c2*Wn+c3*Wn^2;
C4 = c0*W*d/dt*d/dt+c0*Wn/dt^2+c8;

D1 = B2-B1*A2/A1;
D2 = B3-B1*A3/A1;
D3 = B4-B1*A4/A1;

E1 = C2-C1*A2/A1;
E2 = C3-C1*A3/A1;
E3 = C4-C1*A4/A1;
E4 = E1*D3/D1;
E5 = E2- E1*D2/D1;