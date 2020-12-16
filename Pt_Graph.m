pm = 3447E3;
tp =  0.1;
alpha = 2;
t = 0: 0.001 : 0.5
Pt = zeros(1,numel(t));
for i = 1:numel(t)-1
Pt(:,i) = pm*(1-t(i)/tp)*exp(-alpha*t(i)/tp)
plot(t, Pt)
end


