program main
	
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  REAL(8) abd(6,6),M
	dimension u(50000),v(50000),w(50000)

  open(1, file = 'data.csv', status = 'new')

  A11 = 47893916
  A12 = 5268330
  A16 = -0.7287*10.e-9
  A22 = 47893916
  A26 = 0.7287 *10.e-9
  A66 = 7428400

  B11 = 0.7276*10.e-11;
  B12 = 0.7958*10.e-12; 
  B16 = -0.1515*10.e-27; 
  B22 = 0.7276*10.e-11; 
  B26 = 0.1515*10.e-27;
  B66 = 0.7958*10.e-12; 
  
  D11 = 15.332
  D12 = 1.6865
  D16 = -0.1571*10.e-15
  D22 = 15.332
  D26 = 0.1571*10.e-15
  D66 = 2.3780

  write(*,*)'A11= ',A11
  write(*,*)'A12= ',A12
  write(*,*)'B11= ',B11
  write(*,*)'D11= ',D11
  write(*,*)'D12= ',D12


  alf=0.35

  tp=0.0018

  tt=0.1
  dt=0.000002


  
  Pi = 4*ATAN(1.) 

  NM=tt/dt+1
  

  U11 = 1.
  V11 = 1.
  W11 = 1.
   
  a = 0.22
  b = 0.22
  M = 1824*1.e-3
  M = 3.528




  a0=(a*b*M*U11**2)/4.
  b0=(a*b*M*V11**2)/4.
  c0=(a*b*M*W11**2)/4.

  a1=((a**2*A66 + A11*b**2)*Pi**2*U11**2)/(4.*a*b)/a0
  a2=((a**2*A26 + A16*b**2)*Pi**2*U11*V11)/(4.*a*b)/a0
  a3=0./a0
  a4=0./a0
  a5=0./a0

  b1=((a**2*A22 + A66*b**2)*Pi**2*V11**2)/(4.*a*b)/b0
  b2=((a**2*A26 + A16*b**2)*Pi**2*U11*V11)/(4.*a*b)/b0
  b3=0./b0
  b4=0./b0
  b5=0./b0

  c1=(((b**4*D11+a**4*D22+2*a**2*b**2*(D12+2*D66))*Pi**4*W11**2)/(4.*a**3*b**3))/c0
  c2=(8*(B12 - B66)*Pi**2*W11**3)/(3.*a*b)/c0
  c3=((9*a**4*A22 + 2*a**2*(3*A12 + 4*A66)*b**2+9*A11*b**4)*Pi**4*W11**4)/(128.*a**3*b**3)/c0
  c4=0./c0
  c5=0./c0
  c6=0./c0
  c7=0./c0
  c8=((4*a*b*W11)/Pi**2)*25/c0

  un=0.
  vn=0.
  wn=0.

  u(1)=0.
  v(1)=0.
  w(1)=0.

  write(1,10)un,vn,wn,wn,wn


  dt1=1./dt
  dt2=1./dt**2

  t=0.
  do n=1,NM-1

  t=t+dt
  pt=1.0
  pt=(1.-t/tp)*exp(-alf*t/tp)

  if(t.le.0.001) pt = 29000.*t
  if(t.gt.0.001) pt = 29000.- 29000.*(t-0.001)
  if(t.gt.0.002) pt = 0.
	  
  if(t.gt.tp)pt=0.


  www=w(n)

25	AA1=a1+dt2+a0*dt1
  AA2=a2
  AA3=a3+a4*www
  AA4=dt1*un+(dt2+a0*dt1)*u(n)-a5

  BB1=b2
  BB2=dt2+b1+b0*dt1
  BB3=b3+b4*www
  BB4=dt1*vn+(dt2+b0*dt1)*v(n)-b5

  CC1=c4+c6*www
  CC2=c5+c7*www
  CC3=dt2+c1+c2*www+c3*www**2+c0*dt1
  CC4=dt1*wn+(dt2+c0*dt1)*w(n)+c8*pt

  D1=BB2-BB1*AA2/AA1
  D2=BB3-BB1*AA3/AA1
  D3=BB4-BB1*AA4/AA1

  E1=CC2-CC1*AA2/AA1
  E2=CC3-CC1*AA3/AA1
  E3=CC4-CC1*AA4/AA1
  E4=E1*D3/D1
  E5=E2-E1*D2/D1

  w(n+1)=(E3-E4)/E5
  v(n+1)=(D3-D2*w(n+1))/D1
  u(n+1)=(AA4-AA2*v(n+1)-AA3*w(n+1))/AA1

  eps=abs(w(n+1)-www)

  www=w(n+1)

  if(eps.gt.1.0e-12)goto 25

  un=(u(n+1)-u(n))/dt
  vn=(v(n+1)-v(n))/dt
  wn=(w(n+1)-w(n))/dt



  epsx =  0.199839*w(n+1) 

  write(1,10)wn

  enddo

 

10	 format(5e14.7)

  return
  end