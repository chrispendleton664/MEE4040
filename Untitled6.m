F1 = a0*U2 + a1*U0 + a2*V0 + a3*W0 + a4*W0^2 + a5 == 0;
F2 = b0*V2 + b1*V0 + b2*U0 + b3*W0 + b4*W0^2 + b5 == 0;
F3 = c0*W2 + c1*W0 + c2*W0^2 + c3*W0^3 + c4*U0 + c5*V0 + c6*U0*W0 + c8 == 0;
sol = solve(F1, F2, F3);
sol.W