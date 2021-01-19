program main
	
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  REAL(8) abd(6,6),M,d,z(4,4),st(4),sx(4)
  REAL(8) u(500000),v(500000),wbb(500000),wss(500000)
  integer indx(4)

  nb=4
  
  open(1, file = 'data.csv', status = 'new')
  
  A11 = 47893916
  A12 = 5268330
  A16 = 0
  A22 = 47893916
  A26 = 0
  A66 = 7428400
  
  B11 = 0;
  B12 = 0; 
  B16 = 0; 
  B22 = 0; 
  B26 = 0;
  B66 = 0; 
  
  D11 = 15.332
  D12 = 1.6865
  D16 = 0
  D22 = 15.332
  D26 = 0
  D66 = 2.3780		
  
  
  write(*,*)'A11= ',A11
  write(*,*)'A12= ',A12
  write(*,*)'B11= ',B11
  write(*,*)'D11= ',D11
  write(*,*)'D12= ',D12
  write(*,*)'A44= ',A44
  write(*,*)'A55= ',A55
  
  
  alf=0.35
  
  tp=0.0018
  
  tt=0.025
  dt=0.0000021
  
  !C	DOGAL FREKANS VE SONUM ORANLARI GIRILECEK
  
  Pi = 4*ATAN(1.) 
  FR1 = 6800.0 
  FR2 = 6800.0
  FR3 = 673.00
  SORAN = 0.0
  
  NM=tt/dt+1
  
  !C	MALZEME OZELLIKLERI
  U11 = 1.
  V11 = 1.
  WB1 = 1.
  WS1 = 1.
 
  a = 225e-3
  b = 225e-3
	!M = 1800*0.00028*7 (uzunluklar metre)
	!M = (1800*0.0015*2+16*0.025)
  M = (1800*0.00023*2+32*5.08e-3)


  write(*,*)'M= ',M

	!DENKLEM KATSAYILARI ( HER SINIR ÞARTI ÝÇÝN FARKLIDIR)
  !BURADAKÝLER 4 KENAR BASÝT MESNET

  a0=(a*b**9*M*U11**2)/1260.
  b0=(a**9*b*M*V11**2)/1260.
  c0=(a*b*M*WB1**2)/4.
  d0=(a*b*M*WS1**2)/4.

  a1=((3*a**2*A66*b**7 + A11*b**9*Pi**2)*U11**2)/(315.*a)
  a2=(9*a**4*(A12 + A66)*b**4*U11*V11)/Pi**6
  a3=(16*b**3*(b**2*B11 + a**2*(B12 + 2*B66))*(-12 + Pi**2)*U11*WB1)&
  /(3.*a**2*Pi**3)
  a4=0
  a5=(b**3*(A11*b**2*(45 + Pi**4) + a**2*(90*A66 &
  - A12*(-45 + Pi**4)))*U11*WB1**2)/(240.*a**2*Pi)
  a6=(b**3*(A11*b**2*(45 + Pi**4) + a**2*(90*A66 &
  - A12*(-45 + Pi**4)))*U11*WS1**2)/(240.*a**2*Pi)
  a7=(b**3*(A11*b**2*(45 + Pi**4) + a**2*(90*A66 &
  - A12*(-45 + Pi**4)))*U11*WB1*WS1)/(120.*a**2*Pi)
  a8=0

  b1=(9*a**4*(A12 + A66)*b**4*U11*V11)/Pi**6
  b2=((3*a**7*A66*b**2 + a**9*A22*Pi**2)*V11**2)/(315.*b)
  b3=(16*a**3*(a**2*B22 + b**2*(B12 + 2*B66))*(-12 + Pi**2)*V11*WB1)&
  /(3.*b**2*Pi**3)
  b4=0
  b5=(a**3*(a**2*A22*(45 + Pi**4) + b**2*(90*A66 &
  - A12*(-45 + Pi**4)))*V11*WB1**2)/(240.*b**2*Pi)
  b6=(a**3*(a**2*A22*(45 + Pi**4) + b**2*(90*A66 &
  - A12*(-45 + Pi**4)))*V11*WS1**2)/(240.*b**2*Pi)
  b7=(a**3*(a**2*A22*(45 + Pi**4) + b**2*(90*A66 &
  - A12*(-45 + Pi**4)))*V11*WB1*WS1)/(120.*b**2*Pi)
  b8=0

  c1=(16*b**3*(4*b**2*B11 + a**2*(B12 + 2*B66))*(-12 + Pi**2)&
  *U11*WB1)/(3.*a**2*Pi**3)
  c2=(16*a**3*(4*a**2*B22 + b**2*(B12 + 2*B66))&
  *(-12 + Pi**2)*V11*WB1)/(3.*b**2*Pi**3)
  c3=((b**4*D11 + a**4*D22 + 2*a**2*b**2*(D12 + 2*D66))&
  *Pi**4*WB1**2)/(4.*a**3*b**3)
  c4=(a*D22*Pi**4*WB1*WS1)/(4.*b**3)
  c5=(40*B12*Pi**2*WB1**2)/(9.*a*b) + (8*B66*Pi**2*WB1**2)/(9.*a*b)&
  -(16*(B12 + 2*B66)*Pi**2*WB1**2)/(9.*a*b)
  c6=(-8*b*B11*Pi**2*WS1**2)/(9.*a**3) + (8*B12*Pi**2*WS1**2)&
  /(9.*a*b) - (8*a*B22*Pi**2*WS1**2)/(9.*b**3)&
  -(8*B66*Pi**2*WS1**2)/(9.*a*b)
  c7=(-8*b*B11*Pi**2*WB1*WS1)/(9.*a**3) + (16*B12*Pi**2*WB1*WS1)&
  /(3.*a*b)-(8*a*B22*Pi**2*WB1*WS1)/(9.*b**3)-(16*(B12+2*B66)&
  *Pi**2*WB1*WS1)/(9.*a*b)
  c8=(3*(A12 + A66)*b**3*U11*WB1)/(4.*Pi)+(A11*(b**5/60.+(3*b**5)&
  /(4.*Pi**4))*Pi**3*U11*WB1)/(2.*a**2)-(A12*(b**5/60.+(3*b**5)&
  /(4.*Pi**4))*Pi**3*U11*WB1)/(2.*b**2)
  c9=(3*(A12 + A66)*b**3*U11*WS1)/(4.*Pi) + (A11*(b**5/60.+(3*b**5)&
  /(4.*Pi**4))*Pi**3*U11*WS1)/(2.*a**2)-(A12*(b**5/60.+(3*b**5)&
  /(4.*Pi**4))*Pi**3*U11*WS1)/(2.*b**2)
  c10=(3*a**3*(A12 + A66)*V11*WB1)/(4.*Pi)-(A12*(a**5/60.+(3*a**5)&
  /(4.*Pi**4))*Pi**3*V11*WB1)/(2.*a**2) + (A22*(a**5/60.+(3*a**5)&
  /(4.*Pi**4))*Pi**3*V11*WB1)/b**2
  c11=(3*a**3*(A12 + A66)*V11*WS1)/(4.*Pi)-(A12*(a**5/60.+(3*a**5)&
  /(4.*Pi**4))*Pi**3*V11*WS1)/(2.*a**2) + (A22*(a**5/60.+(3*a**5)&
  /(4.*Pi**4))*Pi**3*V11*WS1)/b**2
  c12=(9*a*A22*Pi**4*WB1**3)/(128.*b**3) + (3*A12*Pi**4*WB1**3)/&
  (64.*a*b) + (A66*Pi**4*WB1**3)/(16.*a*b)&
  +(9*A11*b*Pi**4*WB1**3)/(128.*a**3)
  c13=(9*a*A22*Pi**4*WS1**3)/(128.*b**3) + (3*A12*Pi**4*WS1**3)&
  /(64.*a*b)+(A66*Pi**4*WS1**3)/(16.*a*b)+(9*A11*b*Pi**4*WS1**3)&
  /(128.*a**3)
  c14=(27*a*A22*Pi**4*WB1**2*WS1)/(128.*b**3)&
  +(9*A12*Pi**4*WB1**2*WS1)/(64.*a*b)+(3*A66*Pi**4*WB1**2*WS1)&
  /(16.*a*b) + (27*A11*b*Pi**4*WB1**2*WS1)/(128.*a**3)
  c15=(27*a*A22*Pi**4*WB1*WS1**2)/(128.*b**3) &
  +(9*A12*Pi**4*WB1*WS1**2)/(64.*a*b) + (3*A66*Pi**4*WB1*WS1**2)&
  /(16.*a*b) + (27*A11*b*Pi**4*WB1*WS1**2)/(128.*a**3)

  d1=0
  d2=0
  d3=0
  d4=-((a**2*A44*K44 + A55*b**2*K55)*Pi**2*WS1**2)/(4.*a*b)
  d5=(-16*b*B11*Pi**2*WB1**2)/(9.*a**3)-(32*B12*Pi**2*WB1**2)&
  /(9.*a*b)-(16*a*B22*Pi**2*WB1**2)/(9.*b**3)&
  -(16*B66*Pi**2*WB1**2)/(9.*a*b)
  d6=0
  d7=(-16*b*B11*Pi**2*WB1*WS1)/(9.*a**3) - (32*B12*Pi**2*WB1*WS1)&
  /(9.*a*b)-(16*a*B22*Pi**2*WB1*WS1)/(9.*b**3)&
  -(16*B66*Pi**2*WB1*WS1)/(9.*a*b)
  d8=(3*A66*b**3*U11*WB1)/(4.*Pi) + (A11*(b**5/60. + (3*b**5)&
  /(4.*Pi**4))*Pi**3*U11*WB1)/(2.*a**2) + (A12*(b**5/60. &
  +(3*b**5)/(4.*Pi**4))*Pi**3*U11*WB1)/(2.*b**2)
  d9=(3*A66*b**3*U11*WS1)/(4.*Pi) + (A11*(b**5/60. + (3*b**5)&
  /(4.*Pi**4))*Pi**3*U11*WS1)/(2.*a**2) + (A12*(b**5/60. &
  +(3*b**5)/(4.*Pi**4))*Pi**3*U11*WS1)/(2.*b**2)
  d10=(3*a**3*A66*V11*WB1)/(4.*Pi) + (A12*(a**5/60. + (3*a**5)&
  /(4.*Pi**4))*Pi**3*V11*WB1)/(2.*a**2) + (A22*(a**5/60.&
  +(3*a**5)/(4.*Pi**4))*Pi**3*V11*WB1)/(2.*b**2)
  d11=(3*a**3*A66*V11*WS1)/(4.*Pi) + (A12*(a**5/60. + (3*a**5)&
  /(4.*Pi**4))*Pi**3*V11*WS1)/(2.*a**2) + (A22*(a**5/60.&
  +(3*a**5)/(4.*Pi**4))*Pi**3*V11*WS1)/(2.*b**2)
  d12=(-3*a*A22*Pi**4*WB1**3)/(128.*b**3) - (3*A12*Pi**4*WB1**3)&
  /(64.*a*b)+(A66*Pi**4*WB1**3)/(32.*a*b)-(3*A11*b*Pi**4*WB1**3)&
  /(128.*a**3)
  d13=(-3*a*A22*Pi**4*WS1**3)/(128.*b**3) - (3*A12*Pi**4*WS1**3)&
  /(64.*a*b)+(A66*Pi**4*WS1**3)/(32.*a*b)-(3*A11*b*Pi**4*WS1**3)&
  /(128.*a**3)
  d14=(-9*a*A22*Pi**4*WB1**2*WS1)/(128.*b**3)&
  -(9*A12*Pi**4*WB1**2*WS1)/(64.*a*b)+(3*A66*Pi**4*WB1**2*WS1)&
  /(32.*a*b)-(9*A11*b*Pi**4*WB1**2*WS1)/(128.*a**3)
  d15=(-9*a*A22*Pi**4*WB1*WS1**2)/(128.*b**3)&
  -(9*A12*Pi**4*WB1*WS1**2)/(64.*a*b)+(3*A66*Pi**4*WB1*WS1**2)&
  /(32.*a*b) - (9*A11*b*Pi**4*WB1*WS1**2)/(128.*a**3)



!	SONUM KATSAYILARI 

  aS=a0*2*SORAN*FR1
  bS=b0*2*SORAN*FR2
  cS=c0*2*SORAN*FR3
  dS=d0*2*SORAN*FR3

!	KONTROL ICIN YAZDIRMISTIK 

  write(*,*) ' Pi   = ',Pi
  write(*,*) ' a    = ',a 
  write(*,*) ' c1   = ',c1
  write(*,*) ' b1   = ',b1
  
  write(*,*) ' a0   = ',a0
  write(*,*) ' b0   = ',b0
  write(*,*) ' a3   = ',a3
  write(*,*) ' a5   = ',a5
  write(*,*) ' U11  = ',U11


  un=0.
  vn=0.
  wbbn=0.
  wssn=0.
  
  u(1)=0.
  v(1)=0.
  wbb(1)=0.
  wss(1)=0.
  
  write(1,10)un,vn,vn
  
  
  
  dt1=1./dt
  dt2=1./dt**2
  !cccccccccccccccccccccccccc
  t=0.
  do n=1,NM-1
  
  t=t+dt
  !c----------- Renfu LI BLAST LOAD ------------------
  pt=60.86e6*exp(-t/3.33435)    
  !c------------------------------------------
  
  !c-----------  BLAST (EXPONENTIAL) LOAD ------------------
  !c	pt=30e3*(1.-t/tp)*exp(-alf*t/tp)    
  !c------------------------------------------
  
  !c-----------  STEP  LOAD ------------------
  !c	if(t.LE.1*tp) pt=29e3               
  !c	if(t.GT.1*tp) pt=0		             
  !c------------------------------------------
  
  !c-----------  TRIANGULAR LOAD ------------------
  !c	if(t.LE.1*tp) pt=29e3*(1.-t/tp)     
  !c	if(t.GT.1*tp) pt=0					 
  !c------------------------------------------
  
  !c-----------  SINE (HARMONIC) LOAD ---------------
  !c	if(t.LE.1*tp) pt=29e3*SIN(Pi*t/(1*tp))    
  !c	if(t.GT.1*tp) pt=0					 
  !c------------------------------------------
  
  !c-----------  NUCLEAR BLAST (STEPPED TRIANGULAR) LOAD Vinson pg.169-170---------------
  !if(t.LE.0.0009) pt=29e3*(1.-t/tp)     
  !if(t.GT.0.0009) pt=20e3*(1.-t/0.003)	    
  !if(t.GT.0.003) pt=0					 
  !c------------------------------------------
  
  
  c16=(-4*a*b*WB1*Pt)/Pi**2
  d16=(-4*a*b*WS1*Pt)/Pi**2
  
  !CCCCCCCCCCCC
  wb=wbb(n)
  ws=wss(n)
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !cccccccccc  Katsayýlar Matrisi Z cccccccccccccccccccccccccccccccccccccc
  
  25	z(1,1)=(dt2*a0+dt1*aS+a1)
  z(1,2)=(a2)
  z(1,3)=(a3+a5*wb)
  z(1,4)=(a4+a6*ws+a7*wb)
  z(2,1)=(b1)
  z(2,2)=(dt2*b0+dt1*bS+b2)
  z(2,3)=(b3+b5*wb)
  z(2,4)=(b4+b6*ws+b7*wb)
  z(3,1)=(c1+c8*wb+c9*ws)
  z(3,2)=(c2+c10*wb+c11*ws)
  z(3,3)=(dt2*c0+dt1*cS+c3+c5*wb+c7*ws+c12*wb**2+c14*wb*ws)
  z(3,4)=(c4+c6*ws+c13*ws**2+c15*wb*ws)
  z(4,1)=(d1+d8*wb+d9*ws)
  z(4,2)=(d2+d10*wb+d11*ws)
  z(4,3)=(d3+d5*wb+d7*ws+d12*wb**2+d14*wb*ws)
  z(4,4)=(dt2*d0+dt1*dS+d4+d6*ws+d13*ws**2+d15*wb*ws)
  
  ST(1)=dt1*a0*un+(dt2*a0+dt1*aS)*u(n)-a8
  ST(2)=dt1*b0*vn+(dt2*b0+dt1*bS)*v(n)-b8
  ST(3)=dt1*c0*wbbn+(dt2*c0+dt1*cS)*wbb(n)-c16
  ST(4)=dt1*d0*wssn+(dt2*d0+dt1*dS)*wss(n)-d16
  
  !c
  call ludcmp(z,nb,nb,indx,d)
  call lubksb(z,nb,nb,indx,st)
  !c	call gauss(z,st,st,nb)
  
  u(n+1)=st(1)
  v(n+1)=st(2)
  wbb(n+1)=st(3)
  wss(n+1)=st(4)
  
  !c	write(*,*)n,st(1),st(2),st(3),st(4)
  
  eps=abs(wbb(n+1)-wb)+abs(wss(n+1)-ws)
  
  wb=wbb(n+1)
  ws=wss(n+1)
  !c	write(*,*)n,eps
  if(eps.gt.1.0e-40)goto 25
  
  un=(u(n+1)-u(n))/dt
  vn=(v(n+1)-v(n))/dt
  wbbn=(wbb(n+1)-wbb(n))/dt
  wssn=(wss(n+1)-wss(n))/dt
  
  !C	UZAMA HER SINIR SARTI ICIN YAKLASIM FONKSIYONLARINA GORE YENIDEN  HESAPLANMALIDIR 
  !C     BURDAKÝLER 4 KENAR BASIT MESNET

!c	TOP	    a/b= 1.0	X=A/2 Y=B/2        SSSS-New
  epsx = 31.0280755910103*(3.141592653589793*(0.00554*wbb(n+1)))

!c	BOTTOM	a/b= 1.0	X=A/2 Y=B/2        SSSS-New
!c	epsx = 31.0280755910103*(3.141592653589793*(-0.00554*wbb(n+1)))

!c  	write(1,10)t*1000,u(n+1)*1e3,v(n+1)*1e3,
!c	1      w(n+1)*1e3,epsx*1e6
if(mod(n+1,5).eq.0)write(1,10)(wbb(n+1)*1e3)
!c     1  epsx*1e6

!c 	 write(*,10)t*1000,u(n+1)*1e3,v(n+1)*1e3,
!c	1                           wbb(n+1)*1e3,wss(n+1)*1e3
enddo




10	 format(7e14.5)

return
end
!cccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccc
  SUBROUTINE lubksb(a,n,np,indx,b)
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  INTEGER n,np,indx(n)
  REAL(8) a(np,np),b(n)
  INTEGER i,ii,j,ll
  REAL sum
  ii=0
  do 12 i=1,n
    ll=indx(i)
    sum=b(ll)
    b(ll)=b(i)
    if (ii.ne.0)then
      do 11 j=ii,i-1
        sum=sum-a(i,j)*b(j)
11        continue
    else if (sum.ne.0.) then
      ii=i
    endif
    b(i)=sum
12    continue
  do 14 i=n,1,-1
    sum=b(i)
    do 13 j=i+1,n
      sum=sum-a(i,j)*b(j)
13      continue
    b(i)=sum/a(i,i)
14    continue
  return
  END

  SUBROUTINE ludcmp(a,n,np,indx,d)
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  INTEGER n,np,indx(n),NMAX
  REAL(8) d,a(np,np),TINY
  PARAMETER (NMAX=500,TINY=1.0e-20)
  INTEGER i,imax,j,k
  REAL aamax,dum,sum,vv(NMAX)
  d=1.
  do 12 i=1,n
    aamax=0.
    do 11 j=1,n
      if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
    if (aamax.eq.0.) pause 'singular matrix in ludcmp'
    vv(i)=1./aamax
12    continue
  do 19 j=1,n
    do 14 i=1,j-1
      sum=a(i,j)
      do 13 k=1,i-1
        sum=sum-a(i,k)*a(k,j)
13        continue
      a(i,j)=sum
14      continue
    aamax=0.
    do 16 i=j,n

      sum=a(i,j)
      do 15 k=1,j-1
        sum=sum-a(i,k)*a(k,j)
15        continue
      a(i,j)=sum
      dum=vv(i)*abs(sum)
      if (dum.ge.aamax) then
        imax=i
        aamax=dum
      endif
16      continue
    if (j.ne.imax)then
      do 17 k=1,n
        dum=a(imax,k)
        a(imax,k)=a(j,k)
        a(j,k)=dum
17        continue
      d=-d
      vv(imax)=vv(j)
    endif
    indx(j)=imax
    if(a(j,j).eq.0.)a(j,j)=TINY
    if(j.ne.n)then
      dum=1./a(j,j)

      do 18 i=j+1,n
        a(i,j)=a(i,j)*dum
18        continue
    endif
19    continue
  return
  END
!cccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccc
!subroutine gauss(a,c,x,n)
!  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!  REAL(8) a(4,4),c(4),x(4)
!  
!  
!  
!  do k=1,n-1
!  do i=k+1,n
!  m=a(i,k)/a(k,k)
!  c(i)=c(i)-m*c(k)
!  
!    do j=1,k
!    a(i,j)=0
!    enddo
!  
!    do j=k+1,n
!    a(i,j)=a(i,j)-m*a(k,j)
!    enddo
!  enddo
!  enddo
!  
!  
!  do k=n,1,-1
!  top=0.
!    do j=k+1,n
!    top=top+a(k,j)*x(j)
!    enddo
!  
!  x(k)=(c(k)-top)/a(k,k)
!  enddo
!
!
!return
!end