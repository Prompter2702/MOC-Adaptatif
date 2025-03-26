      SUBROUTINE CHKCMX(nc,ns,rv,rs)

      REAL(KIND=8) , DIMENSION(nc,nc) :: rv
      REAL(KIND=8) , DIMENSION(nc,ns) :: rs
      INTEGER      , DIMENSION(2)     :: k

      DO i=1,nc
         rv(i,i)=rv(i,i)-1.0D0
      END DO

      ij=IDAMAX(nc*nc,rv,1)
      CALL ARRPOS(ij,2,(/nc,nc/),k)
      errv=rv(k(1),k(2))

      ij=IDAMAX(nc*ns,rs,1)
      CALL ARRPOS(ij,2,(/nc,ns/),k)
      errs=rs(k(1),k(2))

      WRITE(9,'(1p,2e12.0)')errv,errs

      END SUBROUTINE
