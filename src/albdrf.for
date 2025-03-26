      SUBROUTINE ALBDRF(axis,side,bflx,dirc,w,pisn,nb,ns,nd,ndir,noct)

!     Treats albedo condition on one outer boundary.
!     Incoming angular flux "bflx" for a specified side "side" in all
!     directions is set to constant value.
!     "side" takes value of 1 for minimum and 2 for maximum coordinate
!     of boundary plane.
!     "nb", "ns" and "ndir" are also in input and represent respectively
!     number of space moments, number of surfaces on specified side and
!     number of directions per octant.


      INTEGER       :: axis,side
      REAL          :: bflx(nb,ns,2,nd)
      REAL (KIND=8) :: dirc(ndir),w(ndir),pisn

!     Octant lists: octant numbers of outgoing and incoming directions
!     for two sides "boct(octnber,side,orientation,axis)".

      INTEGER , PARAMETER , DIMENSION(4,2,2,3) :: octl = RESHAPE

     &     ((/1,4,5,8, 2,3,6,7, 2,3,6,7, 1,4,5,8,
     &        1,2,5,6, 3,4,7,8, 3,4,7,8, 1,2,5,6,
     &        1,2,3,4, 5,6,7,8, 5,6,7,8, 1,2,3,4/), (/4,2,2,3/))

      INTEGER , PARAMETER :: incf = 1
      INTEGER , PARAMETER :: outf = 2

      INTEGER       :: d,dd,oc,oct,b,s
      REAL (KIND=8) :: aflx

      DO s=1,ns
      DO b=1,nb

!        Compute total outgoing current.

         aflx=0.
         DO oc=1,noct/2
            oct=octl(oc,side,outf,axis)
            dd=(oct-1)*ndir
            DO d=1,ndir
               dd=dd+1
               aflx=aflx+w(d)*dirc(d)*bflx(b,s,side,dd)
            ENDDO
         ENDDO

!        Set incoming flux.

         aflx=2*aflx/pisn
         DO oc=1,noct/2
            oct=octl(oc,side,incf,axis)
            dd=(oct-1)*ndir
            DO d=1,ndir
               dd=dd+1
               bflx(b,s,side,dd)=aflx
            ENDDO
         ENDDO
      ENDDO
      ENDDO

      END
