      SUBROUTINE COF3C1(ox,oy,oz,sigt,ccof,icof,ecof,tcof)

!     Computes collision "ccof", incoming "icof", escape "ecof"
!     and transmission "tcof" coefficients in one angular direction
!     in one mesh cell for the 3D simplified linear characteristics
!     scheme.

!     Input variables are:

!     "ox"    2 mu  / delta x
!     "oy"    2 eta / delta y
!     "oz"    2 ksi / delta z
!     "sigt"  total cross section.


      IMPLICIT REAL(KIND=8) (A-H,O-Z)

      REAL :: ox,oy,oz,sigt,x,y,z
      REAL , DIMENSION(4,4) :: ccof
      REAL , DIMENSION(4,9) :: icof
      REAL , DIMENSION(9,4) :: ecof
      REAL , DIMENSION(9,9) :: tcof

      REAL(KIND=8) :: lmax

      INTEGER , DIMENSION(3) :: dpos = (/1,2,3/)
      INTEGER , DIMENSION(3) :: ipos = (/1,3,2/)

      INTEGER , DIMENSION(4,3) :: vord = RESHAPE ((/
     &  1,4,3,2, 1,2,4,3, 1,2,3,4/), (/4,3/))

      REAL(KIND=8) , PARAMETER , DIMENSION(3,4) :: ierc = RESHAPE ((/

     &  0.5, -1.5, -1.5,
     & -0.166666666666667,  0.5,  0.5,
     & -0.166666666666667,  0.5,  0.5,
     & -0.166666666666667,  0.5,  0.5

     &  /), (/3,4/))

      REAL , DIMENSION(3) :: xyz
      INTEGER :: b,c,k

      INTEGER , PARAMETER :: nfun=5
      REAL(KIND=8)        :: f(0:nfun)


!     Max. trajectory length, functions.

      lmax=MIN(2./ox,2./oy,2./oz)
         
      x=ox*lmax
      y=oy*lmax
      z=oz*lmax
      t=sigt*lmax

      CALL INCEXP(e,f,t,nfun)


!     Transmission

      CALL CT3C1D(tcof,9,y,z,e,0,0)
      CALL CT3C1O(tcof,9,x,y,z,f,0,3,dpos,ipos)
      CALL CT3C1O(tcof,9,x,z,y,f,0,6,ipos,ipos)

      CALL CT3C1O(tcof,9,y,x,z,f,3,0,dpos,ipos)
      CALL CT3C1D(tcof,9,x,z,e,3,3)
      CALL CT3C1O(tcof,9,y,z,x,f,3,6,ipos,dpos)

      CALL CT3C1O(tcof,9,z,x,y,f,6,0,dpos,dpos)
      CALL CT3C1O(tcof,9,z,y,x,f,6,3,ipos,dpos)
      CALL CT3C1D(tcof,9,x,y,e,6,6)


!     Escape

      CALL CE3C1Z(ecof,9,x,z,y,lmax,f,0,ipos,vord(1,1))
      CALL CE3C1Z(ecof,9,y,x,z,lmax,f,3,dpos,vord(1,2))
      CALL CE3C1Z(ecof,9,z,x,y,lmax,f,6,dpos,vord(1,3))


!     Surface-to-volume (reciprocity)

      xyz(1)=ox
      xyz(2)=oy
      xyz(3)=oz
      DO k=1,3
         DO b=3*(k-1)+1,3*k
         DO c=1,4
            icof(c,b)=xyz(k)*ierc(b-3*k+3,vord(c,k))*ecof(b,c)
         END DO
         END DO
      END DO


!     Collision

      CALL CC3C1A(ccof,4,x,y,z,lmax,f)


!     Conservation check

      CALL CHKC31(lmax,x,y,z,t,ccof,icof,ecof,tcof)

      END
