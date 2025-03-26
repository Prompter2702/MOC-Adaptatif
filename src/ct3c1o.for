      SUBROUTINE CT3C1O(t,nb,z,x,y,f,i,j,ii,jj)

!     Computes a part of transmission matrix "t(nb,nb)" corresponding
!     to an off diagonal block x->z, x.ne.z in simplified linear
!     3D characteristics scheme.
!     Input variables are:
!     z    = oz*lmax (exiting)
!     x    = ox*lmax (entering)
!     y    = oy*lmax (other)
!     f(n) = (1 - Sum[tmax^i/i!,{i,0,n}] Exp[-tmax]) / tmax^(n+1)
!     i,j  = address shift in "t" matrix
!     ii,jj = moments ordering (side dependent)

      REAL , DIMENSION(nb,nb):: t
      REAL(KIND=8) , DIMENSION(0:*)  :: f
      REAL        :: z,x,y
      INTEGER     :: nb,i,j
      INTEGER , DIMENSION(*) :: ii,jj


      t(i+ii(1),j+jj(1))=0.5*x*(f(0)-0.5*y*f(1))

      t(i+ii(2),j+jj(1))=0.5*x*(-f(0)+(x+y/2.)*f(1)-x*y*f(2))

      t(i+ii(3),j+jj(1))=0.25*y*x*(f(1)-y*f(2))

      t(i+ii(1),j+jj(2))=-0.75*y*x*(f(1)-y*f(2))

      t(i+ii(2),j+jj(2))=0.125*3*x*y*(6*x*y*f(3)-2*(2*x+y)*f(2)+2*f(1))

      t(i+ii(3),j+jj(2))=0.125*x*(6*y**3*f(3)-6*y*f(1)+4*f(0))

      t(i+ii(1),j+jj(3))=0.75*x*(2*y*z*f(2)-(y+2*z)*f(1)+2*f(0))

      t(i+ii(2),j+jj(3))=0.75*x*(6*x*y*z*f(3)-2*(y*z+x*(y+2*z))*f(2)
     &                   +(2*x+(y+2*z))*f(1)-2*f(0))

      t(i+ii(3),j+jj(3))=0.375*x*y*(2*f(1)-2*(y+2*z)*f(2)+6*y*z*f(3))

      END SUBROUTINE

