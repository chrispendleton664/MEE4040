function FDECoef = Structure_FDE_coef(NLDE,dt,U0,U1,V0,V1,Wn,W1)
% Coefficients in the finite difference equations
A1 = NLDE.a0/(dt^2) + NLDE.a1;
A2 = NLDE.a2;
A3 = NLDE.a3+ (NLDE.a4*Wn);
A4 = (NLDE.a0*U1)/dt + (NLDE.a0*U0)/(dt^2) - NLDE.a5;

B1 = NLDE.b2;
B2 = (NLDE.b0/dt^2) + NLDE.b1;
B3 = NLDE.b3 + (NLDE.b4*Wn);
B4 = (NLDE.b0*V1)/dt + (NLDE.b0*V0)/(dt^2) - NLDE.b5;

C1 = NLDE.c4 + (NLDE.c6*Wn);
C2 = NLDE.c5 + (NLDE.c7*Wn);
C3 = (NLDE.c0/(dt^2)) + NLDE.c1 + (NLDE.c2*Wn) + NLDE.c3*(Wn^2);
C4 = (NLDE.c0*W1)/dt + (NLDE.c0*Wn)/(dt^2) + NLDE.c8;

D1 = B2 - ((B1*A2)/A1);
D2 = B3 - ((B1*A3)/A1);
D3 = B4 - ((B1*A4)/A1);

E11 = C2 - ((C1*A2)/A1);
E22 = C3 - ((C1*A3)/A1);
E3 = C4 - ((C1*A4)/A1);
E4 = (E11*D3)/D1;
E5 = E22 - ((E11*D2)/D1);
% 
% 	qz = Pt(:,i);
% 
%     dt1=1/dt;
% 	dt2=1/dt^2;
%     
%     
%     AA1=NLDE.a1+dt2+NLDE.a0*dt1;
% 	AA2=NLDE.a2;
% 	AA3=NLDE.a3+NLDE.a4*Wn;
% 	AA4=dt1*U1+(dt2+NLDE.a0*dt1)*Un-NLDE.a5;
% 
% 	BB1=NLDE.b2;
% 	BB2=dt2+NLDE.b1+NLDE.b0*dt1;
% 	BB3=NLDE.b3+NLDE.b4*Wn;
% 	BB4=dt1*V1+(dt2+NLDE.b0*dt1)*Vn-NLDE.b5;
% 
% 	CC1=NLDE.c4+NLDE.c6*Wn;
% 	CC2=NLDE.c5+NLDE.c7*Wn;
% 	CC3=dt2+NLDE.c1+NLDE.c2*Wn+NLDE.c3*Wn^2+NLDE.c0*dt1;
% 	CC4=dt1*W1+(dt2+NLDE.c0*dt1)*Wn+NLDE.c8*qz;
% 
% 	D1=BB2-BB1*AA2/AA1;
% 	D2=BB3-BB1*AA3/AA1;
% 	D3=BB4-BB1*AA4/AA1;
% 
% 	E1=CC2-CC1*AA2/AA1;
% 	E2=CC3-CC1*AA3/AA1;
% 	E3=CC4-CC1*AA4/AA1;
% 	E4=E1*D3/D1;
% 	E5=E2-E1*D2/D1;

FDECoef.A1 = A1;
FDECoef.A2 = A2;
FDECoef.A3 = A3;
FDECoef.A4 = A4;

FDECoef.B1 = B1;
FDECoef.B2 = B2;
FDECoef.B3 = B3;
FDECoef.B4 = B4;

FDECoef.C1 = C1;
FDECoef.C2 = C2;
FDECoef.C3 = C3;
FDECoef.C4 = C4;

FDECoef.D1 = D1;
FDECoef.D2 = D2;
FDECoef.D3 = D3;

FDECoef.E11 = E11;
FDECoef.E22 = E22;
FDECoef.E3 = E3;
FDECoef.E4 = E4;
FDECoef.E5 = E5;
end
