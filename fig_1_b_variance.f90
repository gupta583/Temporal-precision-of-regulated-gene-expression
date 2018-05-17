     PROGRAM pnt
         IMPLICIT NONE
         INTEGER::i,j,k,bin(0:2000),count_f,runs,n,k2,b1,b2,l,ii
         REAL*8 F(20000,2),t0,check,mean,t_max,t1,t2,t,dt,norm
         CHARACTER(50):: flnm
         dt = 0.1
         n = 15                             ! number of molecule
         runs = 10000                        ! number of enemble
         t_max= 30.0
!         t0 = 60.0
 11      bin = 0
         F(:,1) = 0.0 
         F(:,2) = 0.0
!         count_f = n1*runs

        Do ii = 1,runs
        
         IF(ii.lt.10)THEN
           WRITE(flnm,'(a13,I1,a4)'),"../data/file_",ii,".dat"
         ELSE IF(ii.ge.10.and.ii.lt.100) THEN
           WRITE(flnm,'(a13,I2,a4)'),"../data/file_",ii,".dat"
         ELSE IF(ii.ge.100.and.ii.lt.1000)THEN
            WRITE(flnm,'(a13,I3,a4)'),"../data/file_",ii,".dat"
         ELSE IF(ii.ge.1000.and.ii.lt.10000)THEN
            WRITE(flnm,'(a13,I4,a4)'),"../data/file_",ii,".dat"
         ELSE
            WRITE(flnm,'(a13,I5,a4)'),"../data/file_",ii,".dat"
         END IF

         OPEN(11,file = flnm )
         
          i = 0
          t1 = 0.0

          DO while(t1.le.t_max)
            READ(11,*)t1,t2
            i =  i + 1
            F(i,1) = t1 ; F(i,2) = t2
          END DO

       DO j = 1,i
          IF(F(j,2).eq.n)THEN
             b1 = Floor(F(j,1)/dt)
     !        b2 = Floor(F(j+1,1)/dt)
             bin(b1) = bin(b1) + 1
     !        bin(b1:b2) = bin(b1:b2) + 1
             go to 21
     !        IF(b1.eq.0) PRINT*,F(j,1),ii,j
            ! IF(b2.ge.1000) PRINT*,b2,j,f(j,1)
          END IF
       END DO

21       close(11)         
      END Do

         norm = 0.0

        DO i=0,2000
        norm = norm + bin(i)*dt   
        END DO 
        
!        Do i =0,1000,10
!           print*,i*dt,bin(i)/norm
!        END DO

        t2 = 0.0
        mean = 0.0
        DO i =0,2000
            mean = mean + dt*i*dt*bin(i)/norm
            t2   = t2 + ((i*dt)**2)*dt*bin(i)/norm
        END DO

        PRINT*, mean,n,(t2-mean**2)

!        n = n + 1
!        IF(n.le.20) goto 11







      END PROGRAM
