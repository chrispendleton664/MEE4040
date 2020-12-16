% Zafer method
clear all
Q0 = [U0;
      V0;
      W0];
Q2 = [U2;
      V2;
      W2];
M = [a0 0 0;
     0 b0 0;
     0 0 c0];
KL = [a1 a2 a3;
      b2 b1 b3;
      c4 c5 c1];
KNL = [0 0 a4*W0;
       0 0 b4*W0;
       c6*W0 c7*W0 (c2*W0+c3*W0^2)];
F = [a5;
     b5;
     c8];
K = K + KL

     
     