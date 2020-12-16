function Results = NLDE_Solver()
clc;
global Coefficients 
a = 2.54;
b = 2.54;
x = 1.27;
y = 1.27;
h = 0.17;
dt = 0.002E-3;
pm = 3447E3;
tp =  0.1;
alpha = 2;
t = 0:dt:0.006E-2;
Pt = zeros(1,numel(t));
Qq = zeros(3,numel(t));
Qqdot = zeros(3,numel(t));
Qqmax = zeros(3,numel(t));
Qdisp = zeros(3,numel(t));

    for i = 1:numel(t)-1
        Pt(:,i) = pm*(1-t(i)/tp)*exp(-alpha*t(i)/tp);
        Coefficients = GetCoefficients(Qq(:,i), Qqdot(:,i), Pt(:,i));
        fun1 = @root2d;
        x0 = Qq(3,i);
        options = optimoptions('fsolve','FiniteDifferenceType','central','FunctionTolerance',1e-32,'MaxIterations',9999);
        Qq(3,i+1) = fsolve(fun1,x0,options)
%       Qq(3,i+1) = Coefficients.E3-Coefficients.E4/Coefficients.E5;
        Qq(2,i+1) = Coefficients.D3-Coefficients.D2*Qq(3,i)/Coefficients.D1;
        Qq(1,i+1) = 1/Coefficients.A1*(Coefficients.A4-Coefficients.A2*Qq(2,i)-Coefficients.A3*Qq(3,i)); 
%       Qdisp(3,i+1) = lsqr(Qq(3,i+1),Qq(3,i))
%       Qdisp(2,i+1) = lsqr(Qq(2,i+1),Qq(2,i))
%       Qdisp(1,i+1) = lsqr(Qq(1,i+1),Qq(1,i))
        Qqdot(:,i+i) = (Qq(:,i+1) - Qdisp(:,i))/dt;
        Qqmax(1,i+1) = Qq(1,i+1)*sind(2*pi*x/a)*y^2*(y-b)^2;
        Qqmax(2,i+1) = Qq(2,i+1)*sind(2*pi*y/b)*x^2*(x-a)^2;
        Qqmax(3,i+1) = Qq(3,i+1)*sind(pi*x/a)*sin(pi*y/b);
%       wc = Qqmax(:,i)/h;
    end
    plot(t,Qqmax(3,:))
    Results.W = Qqmax(3,:);
%     plot(t, Pt);
%     plot(t, wc);
end

