function Zafer
clc
  A11 = 47893916;
  A12 = 5268330;
  A16 = -0.7287*10.e-9;
  A22 = 47893916;
  A26 = 0.7287*10.e-9;
  A66 = 7428400;

  B11 = 0.7276*10e-11;
  B12 = 0.7958*10e-12; 
  B16 = -0.1515*10e-27; 
  B22 = 0.7276*10e-11; 
  B26 = 0.1515*10e-27;
  B66 = 0.7958*10e-12; 
  
  D11 = 15.332;
  D12 = 1.6865;
  D16 = -0.1571*10.e-15;
  D22 = 15.332;
  D26 = 0.1571*10.e-15;
  D66 = 2.3780;

  alf=0.35;
  tp=0.0018;

  tt=0.1;
  dt=0.000002;
  
  Pi = 180; 

  NM=tt/dt+1;
  NM=int32(NM); 

  U11 = 1;
  V11 = 1;
  W11 = 1;
   
  a = 0.22;
  b = 0.22;
  M = 1824*1.e-3;
  M = 3.528;




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

  un=0;
  vn=0;
  wn=0;
  
  u = zeros(1,(NM-1));
  v = zeros(1,(NM-1));
  w = zeros(1,(NM-1));

  dt1=1/dt;
  dt2=1/dt^2;

  t=0;
  
  for i = 1:(NM)-1

  u(1)=0;
  v(1)=0;
  w(1)=0;    
      
  t=t+dt;
  pt=1.0;
  pt=(1-t/tp)*exp(-alf*t/tp);

  if(t<0.001) pt = 29000*t;
  if(t>0.001) pt = 29000- 29000*(t-0.001);
  if(t>0.002) pt = 0;
  if(t>tp) pt=0;
  end
  end
  end
  end


  www=w(i);
  
  eps = 1;
  
    while (eps > 1.0e-12)
  AA1=a1+dt2+a0*dt1;
  AA2=a2;
  AA3=a3+a4*www;
  AA4=dt1*un+(dt2+a0*dt1)*u(i)-a5;

  BB1=b2;
  BB2=dt2+b1+b0*dt1;
  BB3=b3+b4*www;
  BB4=dt1*vn+(dt2+b0*dt1)*v(i)-b5;

  CC1=c4+c6*www;
  CC2=c5+c7*www;
  CC3=dt2+c1+c2*www+c3*www^2+c0*dt1;
  CC4=dt1*wn+(dt2+c0*dt1)*w(i)+c8*pt;

  D1=BB2-BB1*AA2/AA1;
  D2=BB3-BB1*AA3/AA1;
  D3=BB4-BB1*AA4/AA1;

  E1=CC2-CC1*AA2/AA1;
  E2=CC3-CC1*AA3/AA1;
  E3=CC4-CC1*AA4/AA1;
  E4=E1*D3/D1;
  E5=E2-E1*D2/D1;

  w(i+1)=(E3-E4)/E5;
  v(i+1)=(D3-D2*w(i+1))/D1;
  u(i+1)=(AA4-AA2*v(i+1)-AA3*w(i+1))/AA1;

  eps=abs(w(i+1)-www);

  www=w(i+1);

    end

  un=(u(i+1)-u(i))/dt;
  vn=(v(i+1)-v(i))/dt;
  wn=(w(i+1)-w(i))/dt



  epsx =  0.199839*w(i+1);

  if(mod(i+49,50) == 0)
      fileID = fopen('exp.txt','w');
      fprintf(fileID,'%12s\r\n','w(i+1)');
      fprintf(fileID,'%12.8f\r\n',w(i+1)*1e3);
      fclose(fileID);
%       type exp.txt
  end
% 	if(mod(n+49,50).eq.0)write(1,10)w(n+1)*1e3
%     end
%   
% 
% 
% 10	 format(5e14.7)
  end
%   plot(t, pt);
end  