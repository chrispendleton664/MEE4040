function NLDE = Structure_SSBC_nlde_coef(ABD, Pt, i)

% Laminate Length & Width
% a = 0.22;
% b = 0.22;
% 
% % Blast Pressure Equation
% rho = 1800;
% mdot = rho*a;
% qz = Pt(:,i);
% ABD.B = zeros(3,3);

% a coefficients in the time dependent nonlinear differential equations note:

% a0 = a*b^9*mdot / 1260;
% 
% a1 = (3*a^2*ABD.A(3,3)*b^7+ABD.A(1,1)*b^9*pi^2)/(315*a);
% 
% a2 = 9*a^4*b^4*(ABD.A(1,2)+ABD.A(3,3))/pi^6;
% 
% a3 = 16*b^3*(b^2*ABD.B(1,1)+a^2*(ABD.B(1,2)+2*ABD.B(3,3)))*(-12+pi^2)/(3*a^2*pi^3);
% 
% a4 = b^3*(ABD.A(1,1)*b^2*(45+pi^4)+a^2*(90*ABD.A(3,3)-ABD.A(1,2)*(-45+pi^4)))/(240*a^2*pi);
% 
% a5 = 0;
% 
% % b coefficients in the time dependent nonlinear differential equations 
% 
% b0 = 9*a^9*b*mdot/1260;
% 
% b1 = (3*a^7*b^2*ABD.A(3,3)+a^9*ABD.A(2,2)*pi^2)/(315*b);
% 
% b2 = 9*a^4*b^4*(ABD.A(1,2)+ABD.A(3,3))/pi^6;
% 
% b3 = 16*a^3*(a^2*ABD.B(2,2)+b^2*(ABD.B(1,2)+2*ABD.B(3,3)))*(-12+pi^2)/(3*b^2*pi^3);
% 
% b4 = a^3*(ABD.A(2,2)*a^2*(45+pi^4)+b^2*(90*ABD.A(3,3)-ABD.A(1,2)*(-45+pi^4)))/(240*b^2*pi);
% 
% b5 = 0;
% 
% % c coefficients in the time dependent nonlinear differential equations 
% 
% c0 = a*b*mdot/4;
% 
% c1 = pi^4*(b^4*ABD.D(1,1)+a^4*ABD.D(2,2)+2*a^2*b^2*(ABD.D(1,2)+2*ABD.D(3,3)))/(4*a^3*b^3);
% 
% c2 = 8*pi^2*(ABD.B(1,2)-ABD.B(3,3))/(3*a*b);
% 
% c3 = (9*a^4*ABD.A(2,2)+2*a^2*(3*ABD.A(1,2)+4*ABD.A(3,3))*b^2+9*ABD.A(1,1)*b^4)*pi^4/(128*a^3*b^3);
% 
% c4 = (16*b^3*(4*b^2*ABD.B(1,1)+a^2*(ABD.B(1,2)+2*ABD.B(3,3))))*(-12+pi^2)/(3*a^2*pi^3);
% 
% c5 = (16*a^3*(4*a^2*ABD.B(3,3)+b^2*(ABD.B(1,2)+2*ABD.B(3,3))))*(-12+pi^2)/(3*b^2*pi^3);
% 
% c6 = (3*b^3*(ABD.A(1,2)+ABD.A(3,3)))/(4*pi)  +  (ABD.A(1,1)*((b^5/60)+(3*b^5/(4*pi^4)))*pi^3)/(2*a^2)  -  (ABD.A(1,2)*(b^5/60+3*b^5/(4*pi^4))*pi^3)/(2*b^2); 
% 
% c7 = (3*a^3*(ABD.A(1,2)+ABD.A(3,3)))/(4*pi)  +  (ABD.A(2,2)*((5^5/60)+(3*a^5/(4*pi^4)))*pi^3)/(2*b^2)  -  (ABD.A(1,2)*(a^5/60+3*a^5/(4*pi^4))*pi^3)/(2*a^2);
% 
% c8 = -4*a*b*qz/pi^2;

% A11 = ABD.A(1,1)
% A12 = ABD.A(1,2)
% A16 = ABD.A(1,3)
% A22 = ABD.A(2,2)
% A26 = ABD.A(2,3)
% A66 = ABD.A(3,3)
% 
% B11 = ABD.B(1,1)
% B12 = ABD.B(1,2)
% B16 = ABD.B(1,3)
% B22 = ABD.B(2,2)
% B26 = ABD.B(2,3)
% B66 = ABD.B(3,3)
% 
% D11 = ABD.D(1,1)
% D12 = ABD.D(1,2)
% D16 = ABD.D(1,3)
% D22 = ABD.D(2,2)
% D26 = ABD.D(2,3)
% D66 = ABD.D(3,3)
% 
%       U11 = 1;
%       V11 = 1;
%       W11 = 1;
%      
%     a = 0.5;
% 	b = 0.5;
% 	M = 1824*1.e-3;


  a0=(a*b*M*U11^2)/4;
  b0=(a*b*M*V11^2)/4;
  c0=(a*b*M*W11^2)/4;

  a1=((a^2*A66 + A11*b^2)*Pi^2*U11^2)/(4*a*b)/a0;
  a2=((a^2*A26 + A16*b^2)*Pi^2*U11*V11)/(4*a*b)/a0;
  a3=0/a0;
  a4=0/a0;
  a5=0/a0;

  b1=((a^2*A22 + A66*b^2)*Pi^2*V11^2)/(4*a*b)/b0;
  b2=((a^2*A26 + A16*b^2)*Pi^2*U11*V11)/(4*a*b)/b0;
  b3=0/b0;
  b4=0/b0;
  b5=0/b0;

  c1=(((b^4*D11+a^4*D22+2*a^2*b^2*(D12+2*D66))*Pi^4*W11^2)/(4*a^3*b^3))/c0;
  c2=(8*(B12 - B66)*Pi^2*W11^3)/(3*a*b)/c0;
  c3=((9*a^4*A22 + 2*a^2*(3*A12 + 4*A66)*b^2+9*A11*b^4)*Pi^4*W11^4)/(128*a^3*b^3)/c0;
  c4=0/c0;
  c5=0/c0;
  c6=0/c0;
  c7=0/c0;
  c8=((4*a*b*W11)/Pi^2)*25/c0;




% NLDE Coefficients Structure
NLDE.a0 = a0;
NLDE.a1 = a1;
NLDE.a2 = a2;
NLDE.a3 = a3;
NLDE.a4 = a4;
NLDE.a5 = a5;

NLDE.b0 = b0;
NLDE.b1 = b1;
NLDE.b2 = b2;
NLDE.b3 = b3;
NLDE.b4 = b4;
NLDE.b5 = b5;

NLDE.c0 = c0;
NLDE.c1 = c1;
NLDE.c2 = c2;
NLDE.c3 = c3;
NLDE.c4 = c4;
NLDE.c5 = c5;
NLDE.c6 = c6;
NLDE.c7 = c7;
NLDE.c8 = c8;
end