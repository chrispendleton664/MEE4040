clear all;
E1 = 181E9;
E2 = 10.3E9;
v12 = 0.28;
v21 = 0.28;
G12 = 7.17E9;
fx = 2E6;
fy = -3E6;
fxy = 4E6;
theta1 = 0;
m1 = cosd(theta1);
n1 = sind(theta1);
theta2 = 30;
m2 = cosd(theta1);
n2 = sind(theta1);
theta3 = -45;
m3 = cosd(theta1);
n3 = sind(theta1);
Q11 = E1/(1-(v12*v21))
Q12 = (v21*E1)/(1-(v12*v21))
Q21 = (v12*E2)/((1-v12*v21))
Q22 = E2/(1-(v12*v21))
Q33 = G12
Nx = 1000;
Ny = 1000;
Nxy = 0;
Mx = 0;
My = 0;
Mxy = 0;
n = 1;
z = 0.005
syms p
Qbar111 = (m1^4 * Q11) + (m1^2 * n1^2 * (2*Q12 + 4*Q33)) + (n1^4 * Q22)
Qbar121 =  m1^2 * n1^2 * (Q11 + Q12 - 4*Q33) + Q22*(m1^4 + n1^4)
Qbar131 =  m1^3 * n1 * (Q11 - Q12 - 2*Q33) + m1*n1^3 * (Q12 - Q22 + 2*Q33)
Qbar221 =  n1^4 * Q11 + m1^2 * n1^2*(2*Q12 + 4*Q33) + (m1^4 * Q33)
Qbar231 =  m1 * n1^3 * (Q11 - Q12 - 2*Q33) + m1^3 * n1 *(Q12 - Q22 + 2*Q33)
Qbar331 =  m1^2 * n1^2 * (Q11 + Q12 - 2*Q12 - 2*Q33) + Q33*(m1^4 + n1^4)
QBarMatrix1 = [Qbar111 Qbar121 Qbar131;
              Qbar121 Qbar221 Qbar231;
              Qbar131 Qbar231 Qbar331];
Qbar112 = (m2^4 * Q11) + (m2^2 * n2^2 * (2*Q12 + 4*Q33)) + (n2^4 * Q22)
Qbar122 =  m2^2 * n2^2 * (Q11 + Q12 - 4*Q33) + Q22*(m2^4 + n2^4)
Qbar132 =  m2^3 * n2 * (Q11 - Q12 - 2*Q33) + m2*n2^3 * (Q12 - Q22 + 2*Q33)
Qbar222 =  n2^4 * Q11 + m2^2 * n2^2*(2*Q12 + 4*Q33) + (m2^4 * Q33)
Qbar232 =  m2 * n2^3 * (Q11 - Q12 - 2*Q33) + m2^3 * n2 *(Q12 - Q22 + 2*Q33)
Qbar332 =  m2^2 * n2^2 * (Q11 + Q12 - 2*Q12 - 2*Q33) + Q33*(m2^4 + n2^4)
QBarMatrix2 = [Qbar112 Qbar122 Qbar132;
              Qbar122 Qbar222 Qbar232;
              Qbar132 Qbar232 Qbar332];
Qbar113 = (m3^4 * Q11) + (m3^2 * n3^2 * (2*Q12 + 4*Q33)) + (n3^4 * Q22)
Qbar123 =  m3^2 * n3^2 * (Q11 + Q12 - 4*Q33) + Q22*(m3^4 + n3^4)
Qbar133 =  m3^3 * n3 * (Q11 - Q12 - 2*Q33) + m3*n3^3 * (Q12 - Q22 + 2*Q33)
Qbar223 =  n3^4 * Q11 + m3^2 * n3^2*(2*Q12 + 4*Q33) + (m3^4 * Q33)
Qbar233 =  m3 * n3^3 * (Q11 - Q12 - 2*Q33) + m3^3 * n3 *(Q12 - Q22 + 2*Q33)
Qbar333 =  m3^2 * n3^2 * (Q11 + Q12 - 2*Q12 - 2*Q33) + Q33*(m3^4 + n3^4)
QBarMatrix3 = [Qbar113 Qbar123 Qbar133;
              Qbar123 Qbar223 Qbar233;
              Qbar133 Qbar232 Qbar333];
QBarMatrixTotal = QBarMatrix1+QBarMatrix2+QBarMatrix3
FMatrix = [fx; fy; fxy];
GlobalStrains = linsolve(QBarMatrixTotal,FMatrix)
Aij = symsum(z*p - (z*(p-1)),p,1,n1) * (QBarMatrixTotal)
Bij = symsum(-0.5 * (z * p^2 - z * (p-1)^2),p,1,n1) * (QBarMatrixTotal)
Dij = symsum(1/3 * (z * p^2 - z * (p-1)^2),p,1,n1) * (QBarMatrixTotal)
ABDMatrix = [Aij Bij;
             Bij Dij];
FMMatrix = [Nx;
            Ny;
            Nxy;
            Mx;
            My;
            Mxy];
StrainCurvatures = linsolve(inv(ABDMatrix),FMMatrix)