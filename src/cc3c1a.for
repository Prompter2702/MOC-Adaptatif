      SUBROUTINE CC3C1A(ccof,nc,x,y,z,lmax,f)

!     Computes collision matrix "ccof(nc,nc)" for incoplete linear
!     3D characteristics scheme. Brutal computation without taking
!     into account of reciprocity.
!     Input variables are:
!     x    = ox*lmax
!     y    = oy*lmax
!     z    = oz*lmax

      REAL , DIMENSION(nc,nc) :: ccof
      REAL                    :: x,y,z
      INTEGER                 :: nc

      REAL(KIND=8)                  :: lmax
      REAL(KIND=8) , DIMENSION(0:*) :: f


      ccof(1,1)=0.25*lmax*
     &        (4*f(0)-2*(x+y+z)*f(1)+2*(y*z+x*(y+z))*f(2)-3*x*y*z*f(3))

      ccof(1,2)=-0.25*lmax*3*x*
     &        (2*f(1)-2*(x+y+z)*f(2)+3*(y*z+x*(y+z))*f(3)-6*x*y*z*f(4))

      ccof(1,3)=-0.25*lmax*3*y*
     &        (2*f(1)-2*(x+y+z)*f(2)+3*(y*z+x*(y+z))*f(3)-6*x*y*z*f(4))

      ccof(1,4)=-0.25*lmax*3*z*
     &        (2*f(1)-2*(x+y+z)*f(2)+3*(y*z+x*(y+z))*f(3)-6*x*y*z*f(4))

      ccof(2,1)=0.25*lmax*x*
     &        (2*f(1)-2*(x+y+z)*f(2)+3*(y*z+x*(y+z))*f(3)-6*x*y*z*f(4))

      ccof(2,2)=0.25*lmax*
     &        (4*f(0)-2*(3*x+y+z)*f(1)+(2*y*z+6*x*(y+z))*f(2)+
     &        (6*x**3-9*x*y*z)*f(3)-12*x**3*(y+z)*f(4)+30*x**3*y*z*f(5))

      ccof(2,3)=-0.25*lmax*3*x*y*
     &        (2*f(2)-3*(x+y+z)*f(3)+6*(y*z+x*(y+z))*f(4)-15*x*y*z*f(5))

      ccof(2,4)=-0.25*lmax*3*x*z*
     &        (2*f(2)-3*(x+y+z)*f(3)+6*(y*z+x*(y+z))*f(4)-15*x*y*z*f(5))

      ccof(3,1)=0.25*lmax*y*
     &        (2*f(1)-2*(x+y+z)*f(2)+3*(y*z+x*(y+z))*f(3)-6*x*y*z*f(4))

      ccof(3,2)=-0.25*lmax*3*x*y*
     &        (2*f(2)-3*(x+y+z)*f(3)+6*(y*z+x*(y+z))*f(4)-15*x*y*z*f(5))

      ccof(3,3)=0.25*lmax*
     &        (4*f(0)-2*(x+3*y+z)*f(1)+(6*y*z+2*x*(3*y+z))*f(2)+
     &        (6*y**3-9*x*y*z)*f(3)-12*y**3*(x+z)*f(4)+30*x*y**3*z*f(5))

      ccof(3,4)=-0.25*lmax*3*y*z*
     &        (2*f(2)-3*(x+y+z)*f(3)+6*(y*z+x*(y+z))*f(4)-15*x*y*z*f(5))

      ccof(4,1)=0.25*lmax*z*
     &        (2*f(1)-2*(x+y+z)*f(2)+3*(y*z+x*(y+z))*f(3)-6*x*y*z*f(4))

      ccof(4,2)=-0.25*lmax*3*x*z*
     &        (2*f(2)-3*(x+y+z)*f(3)+6*(y*z+x*(y+z))*f(4)-15*x*y*z*f(5))

      ccof(4,3)=-0.25*lmax*3*y*z*
     &        (2*f(2)-3*(x+y+z)*f(3)+6*(y*z+x*(y+z))*f(4)-15*x*y*z*f(5))

      ccof(4,4)=0.25*lmax*
     &        (4*f(0)-2*(x+y+3*z)*f(1)+(6*y*z+2*x*(y+3*z))*f(2)+
     &        (-9*x*y*z+6*z**3)*f(3)-12*(x+y)*z**3*f(4)+
     &        30*x*y*z**3*f(5))

      END SUBROUTINE
