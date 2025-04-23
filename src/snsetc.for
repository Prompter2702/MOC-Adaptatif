      SUBROUTINE SNSETC(n,m,ityp,mu,et,xi,wg)

      !USE ERRORS
      USE FLGSNQ

      IMPLICIT REAL(KIND=8) (A-H,O-Z)

      REAL(KIND=8) , DIMENSION(*) :: mu,et,xi,wg

      REAL(KIND=8) , ALLOCATABLE , DIMENSION(:)  :: cdir,weig
      REAL(KIND=8) , PARAMETER  :: pio4 = 0.7853981633974483
      LOGICAL                   :: lgdb,lgpr


      lgdb=ityp.EQ.ChebyshevDoubleLegendre .OR.
     &     ityp.EQ.ChebyshevTimesDoubleLegendre

      lgpr=ityp.EQ.ChebyshevTimesLegendre  .OR.
     &     ityp.EQ.ChebyshevTimesDoubleLegendre
      IF(lgpr .AND. n /= m ) THEN
        npolar = m/2
        
      ELSE
        npolar = n/2
      ENDIF 
      
      ALLOCATE(cdir(npolar),weig(npolar),STAT=ierr)

      CALL SNSET1(cdir,weig,lgdb,npolar,.FALSE.)

      IF(lgpr)THEN
         ii=0
         nn=n/2
      ELSE
         ii=1
         nn=0
      ENDIF           
      i=0
      DO k=1,npolar
         r=SQRT(1.0D0-cdir(k)**2)
         nk=nn+ii*k
         DO l=1,nk
            i=i+1
            wg(i)=weig(k)/nk
            xi(i)=cdir(k)
            mu(i)=r*COS((2*nk+1-2*l)*pio4/nk)
            et(i)=SQRT(1.0D0-mu(i)**2-xi(i)**2)
         ENDDO
      ENDDO
      
      DEALLOCATE(cdir,weig,STAT=ierr)

      END
