function F = root2dU(Qq,i,Coefficients)

F(2) = 1/Coefficients.A1*(Coefficients.A4-Coefficients.A2*Qq(2,i)-Coefficients.A3*Qq(3,i)) - Qq(3,i);
end
