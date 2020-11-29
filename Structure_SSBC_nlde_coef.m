ssbcnldecoef = struct(a0, a1, a2, a3, a4, a5, a0, b1, b2, b3, b4, b5, c0, c1, c2, c3, c4, c5, c6, c7, c8)

% coefficients in the time dependent nonlinear differential equations note:
% d1, d2 and d3 are viscous damping coeffiients not sure if needed in this
% case or how to define 

a0 = a*b^9*mdot / 1260;

a1 = 3*a^2*Aijsum(6,6)*b^7+Aijsum(1,1)*b^9*pi^2/(315*a);

a2 = 9*a^4*b^4*(Aijsum(1,2)+Aijsum(6,6))/pi^6;

a3 = 16*b^3*(b^2*Bijsum(1,1)+a^2*(Bijsum(1,2)+2*Bijsum(6,6)))*(-12+pi^2)/(3*a^2*pi^3);

a4 = b^3*(Aijsum(1,1)*b^2*(45+pi^4)+a^2*(90*Aijsum(6,6)-Aijsum(1,2)*(-45+pi^4)))/(240*a^2*pi);

a5 = 0;

%b coefficients in the time dependent nonlinear differential equations 

b0 = 9*a^9*b*mdot/1260;

b1 = 3*a^7*b^2*Aijsum(6,6)+a^9*Aijsum(2,2)*pi^2/(315*b);

b2 = 9*a^4*b^4*(Aijsum(1,2)+Aijsum(6,6))/pi^6;

b3 = 16*a^3*(a^2*Bijsum(2,2)+b^2*(Bijsum(1,2)+2*Bijsum(6,6)))*(-12+pi^2)/(3*b^2*pi^3);

b4 = a^3*(Aijsum(2,2)*a^2*(45+pi^4)+b^2*(90*Aijsum(6,6)-Aijsum(1,2)*(-45+pi^4)))/(240*b^2*pi);

b5 = 0;

%c coefficients in the time dependent nonlinear differential equations 

c0 = a*b*mdot/4;

c1 = pi^4*(b^4*Dijsum(1,1)+A^4*Dijsum(2,2)+2*a^2*b^2*(Dijsum(1,2)+2*Dijsum(6,6)))/(4*a^3*b^3);

c2 = 8*pi^2*(Bijsum(1,2)-Bijsum(6,6))/(3*a*b);

c3 = (9*a^4*Aijsum(2,2)+2*a^2*(3*Aijsum(1,2)+4*Aijsum(6,6))*b^2+9*Aijsum(1,1)*b^4)*pi^4/(128*a^3*b^3);    

c4 = (16*b^3*(4*b^2*Bijsum(1,1)+a^2*(Bijsum(1,2)+2*Bijsum(6,6))))*(-12+pi^2)/(3*a^2*pi^3);

c5 = (16*a^3*(4*a^2*Bijsum(2,2)+b^2*(Bijsum(1,2)+2*Bijsum(6,6))))*(-12+pi^2)/(3*b^2*pi^3);

c6 = (3*b^3*(Aijsum(1,2)+Aijsum(6,6)))/(4*pi)  +  (Aijsum(1,1)*((b^5/60)+(3*b^5/(4*pi^4)))*pi^3)/(2*a^2)  -  (Aijsum(1,2)*(b^5/60+3*b^5/(4*pi^4))*pi^3)/(2*b^2);  

c7 = (3*a^3*(Aijsum(1,2)+Aijsum(6,6)))/(4*pi)  +  (Aijsum(2,2)*((5^5/60)+(3*a^5/(4*pi^4)))*pi^3)/(2*b^2)  -  (Aijsum(1,2)*(a^5/60+3*a^5/(4*pi^4))*pi^3)/(2*a^2);

c8 = -4*a*b*qz/pi^2;