function StrainCurvatures = GUIReferencingMain(MatConst)
syms U V W
% Laminate Engineering Material Constants
% Laminate Force and Moment Instensities
E1 = MatConst.E1;
E2 = MatConst.E2;
v12 = MatConst.v12;
v21 = MatConst.v21;
G12 = MatConst.G12;
NumPlies = MatConst.NumPlies;
fx = MatConst.fx;
fy = MatConst.fy;
fxy = MatConst.fxy;
Nx = MatConst.Nx;
Ny = MatConst.Ny;
Nxy = MatConst.Nxy;
Mx = MatConst.Mx;
My = MatConst.My;
Mxy = MatConst.Mxy;

theta1 = 0;
theta2 = 30;
theta3 = -45;
theta = [theta1,theta2,theta3];
a = 1;
b = 10;
z = [-0.0075, -0.0025, 0.0025, 0.0075];
m = cosd(theta);  
n = sind(theta);

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
% NumPlies
Qb = zeros(3,3,NumPlies);
Aij = zeros(3,3,NumPlies);
Bij = zeros(3,3,NumPlies);
Dij = zeros(3,3,NumPlies);
for i = 1:NumPlies
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
Aijsum = sum(Aij,3);
Bijsum = sum(Bij,3);
Dijsum = sum(Dij,3);
% ABD Matrix 
ABDMatrix = [Aijsum Bijsum;
             Bijsum Dijsum];
% Force and Moment Instensities Matrix
FMMatrix = [Nx; Ny; Nxy; Mx; My; Mxy];
% Laminate Strains and Curvatures
StrainCurvatures = linsolve(ABDMatrix, FMMatrix);

% Solving Time Dependent NLDE

U0 = U;
U1 = diff(U0, 1)
U2 = diff(U0, 2)

V0 = V;
V1 = diff(V0, 1)
V2 = diff(V0, 2)

W0 = W;
W1 = diff(W0, 1)
W2 = diff(W0, 2)
end
