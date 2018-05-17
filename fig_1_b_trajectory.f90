      PROGRAM stochastic
       IMPLICIT NONE
       INTEGER seed1,i,j,A(40000),step,counter,ensemble,B,nstar,seed2
       INTEGER x,c,d,r,d0
       REAL*8 a1,a2,a3,a4,a5,k1,k2,k3,k4,k5
       REAL*8 ran2,r1,tau,it,ft,r2,alpha0,time
       REAL*8 del, eta,tstar,alpha,r0

          
          READ(*,*)seed1
          seed2 = seed1 + 9839
          d  = 0
          x  = 0
          alpha = 4.5
          k5 = 0.1391                      ! Repressor decay
          r0 = 0.4174
          r  = 15                       ! Repressor molecule at t=0
          d0 = 1                        ! c + d = d0
          k1 = alpha/d0                      ! D ==> x+D
          k2 = 0.0                      ! x ==> decay
          k4 = 50.0                      ! C ==> D
          k3 = k4/r0                      ! R+D ==> C
!         seed  = -98392898
          ft    = 30.0
          it    = 0.0
          time  = 0.0

!          nstar = 15
!          tstar = 20.0
!          del   = 2*nstar/tstar
!          eta   = 2.0
!          h = 2.0



          WRITE(*,*) time, x,r,d

          DO while(time.lt.ft)
              
               r1  = ran2(seed1)
               r2  = ran2(seed2)
               a1  = k1*d
               a2  = k2*x
               a3  = k3*r*d
               a4  = k4*(d0-d)
               a5  = k5*r

               alpha0 = a1 + a2 + a3 + a4 + a5
               tau = log(1/r1)/alpha0

               IF(r2.lt.a1/alpha0)THEN
                    x = x + 1
               ELSE IF(r2.lt.(a1+a2)/alpha0)THEN 
                    x = x - 1
               ELSE IF(r2.lt.(a1+a2+a3)/alpha0)THEN
                   d  = d - 1 
               ELSE IF(r2.lt.(a1+a2+a3+a4)/alpha0) THEN
                  d  = d + 1
               else 
                  r  = r - 1 
               END IF

               If(x.lt.0.or.r.lt.0.or.d.lt.0)THEN
               PRINT*, x,r,d
               STOP
               end if

               time = time + tau
               WRITE(*,*) time, x,r,d

          END DO
      
      END PROGRAM stochastic 









!------------------------------------------------------------------------------------------------------------------------
!                                  Function to generate a random number
!------------------------------------------------------------------------------------------------------------------------  
      FUNCTION ran2(idum)
      IMPLICIT NONE
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real*8  ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1& 
        ,IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791&
        ,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/

      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue  
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END  FUNCTION ran2
