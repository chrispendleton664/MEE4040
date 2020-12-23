clear all;
E1 = 181E9;
E2 = 10.3E9;
v12 = 0.28;
v21 = (v12*E2)/E1;
G12 = 7.17E9;
fx = 2E6;
fy = -3E6;
fxy = 4E6;
theta = 60;
Q11 = E1/(1-(v12*v21))
Q12 = (v21*E1)/(1-(v12*v21))
Q21 = (v12*E2)/((1-v12*v21))
Q22 = E2/(1-(v12*v21))
Q33 = G12
m = cosd(theta)
n = sind(theta) 
Qbar11 = (m^4 * Q11) + (m^2 * n^2 * (2*Q12 + 4*Q33)) + (n^4 * Q22)
Qbar12 =  m^2 * n^2 * (Q11 + Q12 - 4*Q33) + Q22*(m^4 + n^4)
Qbar13 =  m^3 * n * (Q11 - Q12 - 2*Q33) + m*n^3 * (Q12 - Q22 + 2*Q33)
Qbar22 =  n^4 * Q11 + m^2 * n^2*(2*Q12 + 4*Q33) + (m^4 * Q33)
Qbar23 =  m * n^3 * (Q11 - Q12 - 2*Q33) + m^3 * n *(Q12 - Q22 + 2*Q33)
Qbar33 =  m^2 * n^2 * (Q11 + Q12 - 2*Q12 - 2*Q33) + Q33*(m^4 + n^4)
A = [Qbar11 Qbar12 Qbar13;
     Qbar12 Qbar22 Qbar23;
     Qbar13 Qbar23 Qbar33];
B = [fx; fy; fxy];
X = linsolve(A,B)