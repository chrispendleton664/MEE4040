format shortEng;
E1 = 181E9;
E2 = 10.3E9;
v12 = 0.28;
v21 = 0.28;
G12 = 7.17E9;
fx = 0;
fy = 0;
fxy = 0;
theta1 = 0;
m1 = cosd(theta1);
n1 = sind(theta1);
theta2 = 30;
m2 = cosd(theta2);
n2 = sind(theta2);
theta3 = -45;
m3 = cosd(theta3);
n3 = sind(theta3);
Nx = 1000;
Ny = 1000;
Nxy = 0;
Mx = 0;
My = 0;
Mxy = 0;
n = 2;
z = 0.005
syms p

S11 = 1/E1 
S12 = -v12/E1
S21 = -v12/E1
S22 = 1/E2
S33 = 1/G12
qwerty = [S11 S12 0;
          S12 S22 0;
          0 0 S33];
Qmatrixcheck = inv(qwerty)
S1Transform = [m1^4 n1^4 2*(m1^2)*n1^2 (m1^2)*n1^2;
              n1^4 m1^4 2*(m1^2)*n1^2 (m1^2)*n1^2;
              4*(m1^2)*n1^2 4*(m1^2)*n1^2 -8*(m1^2)*n1^2 ((m1^2)-n1^2)^2;
              (m1^2)*n1^2 (m1^2)*n1^2 (m1^4)+n1^4 -(m1^2)*n1^2;
              2*(m1^3)*n1 -2*m1*(n1^3) 2*(m1*(n1^3)-(m1^3)*n1) m1*(n1^3)-(m1^3)*n1;
              2*m1*(n1^3) -2*n1*(m1^3) 2*(n1*(m1^3)-(n1^3)*m1) n1*(m1^3)-(n1^3)*m1]
S2Transform = [m2^4 n2^4 2*(m2^2)*n2^2 (m2^2)*n2^2;
              n2^4 m2^4 2*(m2^2)*n2^2 (m2^2)*n2^2;
              4*(m2^2)*n2^2 4*(m2^2)*n2^2 -8*(m2^2)*n2^2 ((m2^2)-n2^2)^2;
              (m2^2)*n2^2 (m2^2)*n2^2 (m2^4)+n2^4 -(m2^2)*n2^2;
              2*(m2^3)*n2 -2*m2*(n2^3) 2*(m2*(n2^3)-(m2^3)*n2) m2*(n2^3)-(m2^3)*n2;
              2*m2*(n2^3) -2*n2*(m2^3) 2*(n2*(m2^3)-(n2^3)*m2) n2*(m2^3)-(n2^3)*m2]
S3Transform = [m3^4 n3^4 2*(m3^2)*n3^2 (m3^2)*n3^2;
              n3^4 m3^4 2*(m3^2)*n3^2 (m3^2)*n3^2;
              4*(m3^2)*n3^2 4*(m3^2)*n3^2 -8*(m3^2)*n3^2 ((m3^2)-n3^2)^2;
              (m3^2)*n3^2 (m3^2)*n3^2 (m3^4)+n3^4 -(m3^2)*n3^2;
              2*(m3^3)*n3 -2*m3*(n3^3) 2*(m3*(n3^3)-(m3^3)*n3) m3*(n3^3)-(m3^3)*n3;
              2*m3*(n3^3) -2*n3*(m3^3) 2*(n3*(m3^3)-(n3^3)*m3) n3*(m3^3)-(n3^3)*m3]
SMatrix = [S11;
           S22;
           S12;
           S33]
S1Bar = S1Transform * SMatrix
S1Full = [S1Bar(1,1) S1Bar(4,1) S1Bar(5,1);
          S1Bar(4,1) S1Bar(2,1) S1Bar(6,1);
          S1Bar(5,1) S1Bar(6,1) S1Bar(3,1)]
QBar1 = inv(S1Full)      
S2Bar = S2Transform * SMatrix
S2Full = [S2Bar(1,1) S2Bar(4,1) S2Bar(5,1);
          S2Bar(4,1) S2Bar(2,1) S2Bar(6,1);
          S2Bar(5,1) S2Bar(6,1) S2Bar(3,1)]   
QBar2 = inv(S2Full)      
S3Bar = S3Transform * SMatrix
S3Full = [S3Bar(1,1) S3Bar(4,1) S3Bar(5,1);
          S3Bar(4,1) S3Bar(2,1) S3Bar(6,1);
          S3Bar(5,1) S3Bar(6,1) S3Bar(3,1)] 
QBar3 = inv(S3Full)      
A = QBar1 + QBar2 + QBar3;
B = [fx; fy; fxy];
Strain = linsolve(A,B)
Aij = QBar1*((-0.0025) - (-0.0075)) + (QBar2*(0.0025 - (-0.0025))) + (QBar3*(0.0075 - 0.0025))
Bij = 0.5*QBar1*((-0.0025)^2 - (-0.0075)^2) + 0.5*(QBar2*(0.0025^2 - (-0.0025)^2)) + 0.5*(QBar3*(0.0075^2 - 0.0025^2))
Dij = (1/3)*QBar1*((-0.0025)^3 - (-0.0075)^3) + (1/3)*(QBar2*(0.0025^3 - (-0.0025)^3)) + (1/3)*(QBar3*(0.0075^3 - 0.0025^3))
ABDMatrix = [Aij Bij;
             Bij Dij];
FMMatrix = [Nx; Ny; Nxy; Mx; My; Mxy];
StrainCurvatures = linsolve(ABDMatrix, FMMatrix)