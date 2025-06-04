      SUBROUTINE INNERR(er,it,flxm,flxp,nr,nh,nc,lgpr,lgeh,negf,ndir,
     &                  iopt)

!     Pointwise convergence test on flux moments.

!     er   is returned value of relative error,
!     it   is input as iteration count,
!     flxp is flux difference in current-previous iteration,
!     flxm is flux in current iteration,
!     nr   is number of regions,
!     nh   is number of flux harmonics,
!     nc   is numer of flux spatial moments,
!     lgpr is .TRUE. if printout,
!     negf is number of negative angular fluxes encoutered,
!     ndir is half total number of discrete directions,
!     iopt is the flag for the error norm 

      REAL    :: er
      INTEGER :: it,nr,nh,nc,negf,ndir
      LOGICAL :: lgpr,lgeh
      REAL     , DIMENSION(nr,nh,nc) :: flxm,flxp

      INTEGER :: r,h,c,ir,ih,ic
      REAL    :: e,aux
      
      real(kind=8) :: norm,norx 
      
      
      
      SELECT CASE (iopt)
      
      CASE( 1 )
      ! SCALAR-FLUX ONLY 
      
      er=0
      ir=0
      ih=0
      ic=0
      DO r=1,nr
         aux=ABS(flxm(r,1,1))
         IF(aux.LE.0.0E0)CYCLE 
            e=ABS(flxp(r,1,1))/aux
            IF(e.LE.er)CYCLE
            er=e
            ir=r
            ih=1
            ic=1
      ENDDO

      IF(negf.GT.0)THEN
         negf=50.*negf/REAL(ndir*nr)
         IF(lgpr)WRITE(9,'(8X,I8,1P,E15.4,I12,3I6)')it,er,ir,ih,ic,negf
         IF(lgeh)WRITE(6,'(8X,I8,1P,E15.4,I12,3I6)')it,er,ir,ih,ic,negf
      ELSE
         IF(lgpr)WRITE(9,'(8X,I8,1P,E15.4,I12,2I6)')it,er,ir,ih,ic
         IF(lgeh)WRITE(6,'(8X,I8,1P,E15.4,I12,2I6)')it,er,ir,ih,ic
      ENDIF
      
      RETURN
      
      CASE( 2 )
      
      ! RMS (Absolute)
       
      norm = 0.0
      DO c=1,nc 
      DO h=1,nh
      DO r=1,nr
         norm = norm + (flxp(r,h,c))**2
         norx = norx + (flxm(r,h,c))**2
      ENDDO
      ENDDO
      ENDDO
      IF(norx/=0.0D0)THEN 
         er =  SQRT(norm/norx)
      ELSE 
         er =  SQRT(norm)
      ENDIF 
      IF(lgpr)WRITE(9,'(8X,I8,1P,E15.4)')it,er
      IF(lgeh)WRITE(6,'(8X,I8,1P,E15.4)')it,er
      RETURN 
      
      CASE( 3 )
      
      
      ! AVERAGE
       
      norm = 0.0
      DO c=1,nc 
      DO h=1,nh
      DO r=1,nr
         aux = MAX(flxm(r,h,c),flxm(r,1,1))
         IF(flxm(r,h,c)<=0.0) aux = 1.0
         norm = norm + ABS(flxp(r,h,c))/aux
      ENDDO
      ENDDO
      ENDDO
      
      er =  norm/(nr*nh*nc)
      IF(lgpr)WRITE(9,'(8X,I8,1P,E15.4)')it,er
      IF(lgeh)WRITE(6,'(8X,I8,1P,E15.4)')it,er
      RETURN
      
      CASE( 4 )
      
      
      ! RMS
       
      norm = 0.0
      DO c=1,nc 
      DO h=1,nh
      DO r=1,nr
         aux = MAX(flxm(r,h,c),flxm(r,1,1))
         IF(flxm(r,h,c)<=0.0) aux = 1.0
         norm = norm + (flxp(r,h,c)/aux)**2
      ENDDO
      ENDDO
      ENDDO
      
      er =  SQRT(norm/(nr*nh*nc))
      IF(lgpr)WRITE(9,'(8X,I8,1P,E15.4)')it,er
      IF(lgeh)WRITE(6,'(8X,I8,1P,E15.4)')it,er
      RETURN 
      
      CASE DEFAULT 
      
      er=0
      ir=0
      ih=0
      ic=0
      DO r=1,nr
         aux=ABS(flxm(r,1,1))
         IF(aux.LE.0.0E0)CYCLE 
         DO h=1,nh
         DO c=1,nc
            e=ABS(flxp(r,h,c))/MAX(aux,ABS(flxm(r,h,c)))
            IF(e.LE.er)CYCLE
            er=e
            ir=r
            ih=h
            ic=c
         ENDDO
         ENDDO
      ENDDO

      IF(negf.GT.0)THEN
         negf=50.*negf/REAL(ndir*nr)
         IF(lgpr)WRITE(9,'(8X,I8,1P,E15.4,I12,3I6)')it,er,ir,ih,ic,negf
         IF(lgeh)WRITE(6,'(8X,I8,1P,E15.4,I12,3I6)')it,er,ir,ih,ic,negf
      ELSE
         IF(lgpr)WRITE(9,'(8X,I8,1P,E15.4,I12,2I6)')it,er,ir,ih,ic
         IF(lgeh)WRITE(6,'(8X,I8,1P,E15.4,I12,2I6)')it,er,ir,ih,ic
      ENDIF
      
      RETURN 
      
      END SELECT 

      END

      SUBROUTINE INNERR_A(er,it,flxm,flxp,adx,vox,rex,
     &                    nbox,nr,nh,nc,lgpr,lgeh,negf,ndir)

!     Pointwise convergence test on flux moments.

!     er   is returned value of relative error,
!     it   is input as iteration count,
!     flxp is flux difference in current-previous iteration,
!     flxm is flux in current iteration,
!     nr   is number of regions,
!     nh   is number of flux harmonics,
!     nc   is numer of flux spatial moments,
!     lgpr is .TRUE. if printout,
!     negf is number of negative angular fluxes encoutered,
!     ndir is half total number of discrete directions,
      IMPLICIT NONE 
      REAL    :: er
      INTEGER :: it,nr,nh,nc,negf,ndir,nbox
      LOGICAL :: lgpr,lgeh
      REAL     , DIMENSION(nr,nh,nc) :: flxm,flxp

      INTEGER :: r,h,c,ir,ih,ic,rr,i
      INTEGER :: adx(*),rex(*)
      REAL    :: vox(*)
      REAL    :: e,fp(nh,nc),fm(nh,nc),aux


      er=0
      ir=0
      ih=0
      ic=0
      DO i=1,nbox
       fp = 0.0
       fm = 0.0
       DO rr=adx(i),adx(i+1)-1
         r=rex(rr) 
         DO h=1,nh
         DO c=1,nc
            fp(h,c) = fp(h,c) + flxp(r,h,c) * vox(rr) 
            fm(h,c) = fm(h,c) + flxm(r,h,c) * vox(rr) 
         ENDDO
         ENDDO
       ENDDO
       aux=ABS(fm(1,1))
       IF(aux.EQ.0.)CYCLE
       DO h=1,nh
       DO c=1,nc
          e=ABS(fp(h,c))/MAX(aux,ABS(fm(h,c)))
          IF(e.LE.er)CYCLE
          er=e
          ir=i
          ih=h
          ic=c
       ENDDO
       ENDDO
      ENDDO

      IF(negf.GT.0)THEN
         negf=50.*negf/REAL(ndir*nr)
         IF(lgpr)WRITE(9,'(8X,I8,1P,E15.4,4I6)')it,er,ir,ih,ic,negf
         IF(lgeh)WRITE(6,'(8X,I8,1P,E15.4,4I6)')it,er,ir,ih,ic,negf
      ELSE
         IF(lgpr)WRITE(9,'(8X,I8,1P,E15.4,3I6)')it,er,ir,ih,ic
         IF(lgeh)WRITE(6,'(8X,I8,1P,E15.4,3I6)')it,er,ir,ih,ic
      ENDIF
      END
      
      SUBROUTINE INNERR_T(er,it,flxm,flxp,nr,nh,nc,lgpr,lgeh,negf,ndir)

!     Pointwise convergence test on flux moments.

!     er   is returned value of relative error,
!     it   is input as iteration count,
!     flxp is flux difference in current-previous iteration,
!     flxm is flux in current iteration,
!     nr   is number of regions,
!     nh   is number of flux harmonics,
!     nc   is numer of flux spatial moments,
!     lgpr is .TRUE. if printout,
!     negf is number of negative angular fluxes encoutered,
!     ndir is half total number of discrete directions,

      REAL    :: er
      INTEGER :: it,nr,nh,nc,negf,ndir
      LOGICAL :: lgpr,lgeh
      REAL     , DIMENSION(nr,nh,nc) :: flxm,flxp

      INTEGER :: r,h,c,ir,ih,ic
      REAL(KIND=8)    :: e,n


      er=0
      ir=0
      ih=0
      ic=0
      e=0.0D0
      n=0.0D0
      DO r=1,nr
         DO h=1,nh
         DO c=1,nc
            e=e+(flxp(r,h,c))**2
            n=n+(flxm(r,h,c))**2
         ENDDO
         ENDDO
      ENDDO
      er = e/n
      IF(negf.GT.0)THEN
         negf=50.*negf/REAL(ndir*nr)
         IF(lgpr)WRITE(9,'(8X,I8,1P,E15.4,4I6)')it,er,ir,ih,ic,negf
         IF(lgeh)WRITE(6,'(8X,I8,1P,E15.4,4I6)')it,er,ir,ih,ic,negf
      ELSE
         IF(lgpr)WRITE(9,'(8X,I8,1P,E15.4,3I6)')it,er,ir,ih,ic
         IF(lgeh)WRITE(6,'(8X,I8,1P,E15.4,3I6)')it,er,ir,ih,ic
      ENDIF

      END
