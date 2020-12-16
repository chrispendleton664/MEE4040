A = delsq(numgrid('L',400));
b = ones(size(A,1),1);
x = pcg(A,b,[],1000);
norm(b-A*x)
Qq = zeros(3,numel(t));
Qq1 = zeros(3,numel(t));
Qq2 = zeros(3,numel(t));

while tol>.00001  % *desired accuracy*
                c = c+1
                Qq1(3,c)=Coefficients.c2*Qq1(3,c)-Coefficients.c3*Qq1(3,c)  % *added*
                Qq(3,c)=Coefficients.E3-Coefficients.E4/Coefficients.E5
                Qq2(3,c)=Coefficients.c2*Qq(3,c)-Coefficients.c3*Qq(3,c)  % * note changes to alpha & c
                tol = abs(Qq2(3,1)-Qq1(3,c))
        end
        while tol>.00001  % *desired accuracy*
                k = k+1;
                Qq(2,k+1)=Coefficients.D3-Coefficients.D2*Qq(3,i)/Coefficients.D1;  % * note changes to alpha & c
                Qq(2,k)=0; % *added*
                tol = abs(Qq(k)-Qq(k+1));
        end
        while tol>.00001  % *desired accuracy*
                c = g+1;
                Qq(1,g+1)=1/Coefficients.A1*(Coefficients.A4-Coefficients.A2*Qq(2,i)-Coefficients.A3*Qq(3,i));  % * note changes to alpha & c
                Qq(1,g)=0; % *added*
                tol = abs(Qq(g)-Qq(g+1));
        end