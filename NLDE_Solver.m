function NLDE_Solver()
clc
a = 0.22;
b = 0.22;
x = 0.11;
y = 0.11;
dt = 0.002;
h = 1.98E-3;
pm = 28.9E3;
tp =  0.0018;
alpha = 0.35;
t = 0:dt:0.03
Pt = zeros(1,numel(t));
Qq = zeros(3,numel(t));
Qqdot = zeros(3,numel(t));
Qqmax = zeros(3,numel(t));
wc = zeros(3,numel(t));
ABD = ABD_and_Strain;

U0 = Qq(1,1);
U1 = Qqdot(1,1);

V0 = Qq(2,1);
V1 = Qqdot(2,1);

Wn = Qq(3,1);
W1 = Qqdot(3,1);

    for i = 1:numel(t)-1
        Pt(:,i) = -pm*(1-t(i)/tp)*exp(-alpha*t(i)/tp);
        NLDE = Structure_SSBC_nlde_coef(ABD, Pt, i);
        diff = 10;
        
        % U, V & W and U, V & W dot
    
            
        while diff > 1E-12 
              
            FDECoef = Structure_FDE_coef(NLDE,dt,U0,U1,V0,V1,Wn,W1);
            
            
            Qq(3,i+1) = (FDECoef.E3-FDECoef.E4)/FDECoef.E5;
            Qq(2,i+1) = (FDECoef.D3-FDECoef.D2*Qq(3,i+1))/FDECoef.D1;
            Qq(1,i+1) = (FDECoef.A4-FDECoef.A2*Qq(2,i+1)-FDECoef.A3*Qq(3,i+1))/FDECoef.A1;

            diff = abs(Qq(3,i+1) - Wn);

            Wn = Qq(3,i+1)
            V0 = Qq(2,i+1)
            U0 = Qq(1,i+1)
            
         end
            
%       Coefficients = GetCoefficients(Qq(:,i), Qqdot(:,i), Pt(:,i));
%         U0 = Qq(1,i);
%         V0 = Qq(2,i);

        Qqdot(:,i) = (Qq(:,i+1) - Qq(:,i))/dt;

%         U1 = Qqdot(1,i);
%         V1 = Qqdot(2,i);
%         W1 = Qqdot(3,i); 
        
        Qqmax(1,i+1) = Qq(1,i+1)*sind(2*pi*x/a)*y^2*(y-b)^2;
        Qqmax(2,i+1) = Qq(2,i+1)*sind(2*pi*y/b)*x^2*(x-a)^2;
        Qqmax(3,i+1) = Qq(3,i+1)*sind(pi*x/a)*sind(pi*y/b);
        wc(:,i) = Qqmax(:,i) / h;
    end
%     plot(t,Qq(3,:))
%     plot(t,Qqmax(3,:))
    plot(t, Pt);
%     plot(t, wc(3,:));
end

