function ABD = ABD_and_Strain
% Laminate Engineering Material Constants
E1 = 181E9;
E2 = 10.3E9;
v12 = 0.28;
v21 = (E2*v12)/E1;
G12 = 7.17E9;
fx = 0;
fy = 0;
fxy = 0;
theta1 = 0;
theta2 = 30;
theta3 = -45;
numPlies = 3;
theta = [theta1,theta2,theta3];
a = 1;
b = 10;
z = [-0.0075, -0.0025, 0.0025, 0.0075];
m = cosd(theta);  
n = sind(theta);


% Laminate Force and Moment Instensities
Nx = 1000;
Ny = 1000;
Nxy = 0;
Mx = 0;
My = 0;
Mxy = 0;

% Q Values
Q11 = E1/(1-(v12*v21));
Q12 = (v21*E1)/(1-(v12*v21));
Q21 = (v12*E2)/((1-v12*v21));
Q22 = E2/(1-(v12*v21));
Q33 = G12;

% Q Matrix
Q =  [Q11 Q12 0;
      Q12 Q22 0;
      0   0   Q33];

% Q Tranformation Matrix, using 3 x 3 x i matrix, with i being 1 to
% numPlies
Qb = zeros(3,3,numPlies);
Aij = zeros(3,3,numPlies);
Bij = zeros(3,3,numPlies);
Dij = zeros(3,3,numPlies);
for i = 1:numPlies
    Qt1 = [m(i)^2 n(i)^2 -2*m(i)*n(i);
          n(i)^2 m(i)^2 2*m(i)*n(i);
          m(i)*n(i) -m(i)*n(i) m(i)^2 - n(i)^2];
    Qt2 = [m(i)^2 n(i)^2 m(i)*n(i);
          n(i)^2 m(i)^2 -m(i)*n(i);
          -2*m(i)*n(i) 2*m(i)*n(i) m(i)^2 - n(i)^2];
    
% Qbar Matrix for each Ply
Qb(:,:,i) = Qt1*Q*Qt2;

% ABD Matricies for each Ply
Aij(:,:,i) = Qb(:,:,i)*(z(i+1)-z(i));
Bij(:,:,i) = Qb(:,:,i)*0.5*(z(i+1)^2 - z(i)^2);
Dij(:,:,i) = Qb(:,:,i)*(1/3)*(z(i+1)^3 - z(i)^3);
end

% Qbar Matricies sum for all plies
AQbarSum = sum(Qb,3);
% Force matrix
B = [fx; fy; fxy];
% In-plane Strains
Strain = linsolve(AQbarSum,B);

% ABD Matricies sum for all plies
ABD.Aijsum = sum(Aij,3)
ABD.Bijsum = sum(Bij,3)
ABD.Dijsum = sum(Dij,3)
% ABD Matrix 
ABDMatrix = [ABD.Aijsum ABD.Bijsum;
             ABD.Bijsum ABD.Dijsum];
% Force and Moment Instensities Matrix
FMMatrix = [Nx; Ny; Nxy; Mx; My; Mxy];
% Laminate Strains and Curvatures
StrainCurvatures = linsolve(ABDMatrix, FMMatrix)
end
