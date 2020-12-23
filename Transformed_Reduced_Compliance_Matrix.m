clear all;
E1 = 181E9;
E2 = 10.3E9;
v12 = 0.28;
v21 = 0.28;
G12 = 7.17E9;
fx = 2E6;
fy = -3E6;
fxy = 4E6;
theta = 60;
m = cosd(theta)
n = sind(theta)
S11 = 1/E1 
S12 = -v12/E1
S21 = -v12/E1
S22 = 1/E2
S33 = 1/G12
STransform = [m^4 n^4 2*(m^2)*n^2 (m^2)*n^2;
              n^4 m^4 2*(m^2)*n^2 (m^2)*n^2;
              4*(m^2)*n^2 4*(m^2)*n^2 -8*(m^2)*n^2 ((m^2)-n^2)^2;
              (m^2)*n^2 (m^2)*n^2 (m^4)+n^4 -(m^2)*n^2;
              2*(m^3)*n -2*m*(n^3) 2*(m*(n^3)-(m^3)*n) m*(n^3)-(m^3)*n;
              2*m*(n^3) -2*n*(m^3) 2*(n*(m^3)-(n^3)*m) n*(m^3)-(n^3)*m]
SMatrix = [S11;
           S22;
           S12;
           S33]
SBar = STransform * SMatrix
SFull = [SBar(1,1) SBar(4,1) SBar(5,1);
         SBar(4,1) SBar(2,1) SBar(6,1);
         SBar(5,1) SBar(6,1) SBar(3,1)]
QBar = inv(SFull)
A = QBar;
B = [fx; fy; fxy];
X = linsolve(A,B)