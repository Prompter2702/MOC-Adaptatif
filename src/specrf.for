      SUBROUTINE SPECRF(bflx,nb,ns,ndir,side,oct,rdir)

!     Treats specular reflection on one outer boundary.
!     Incoming angular flux "bflx" for a specified side "side" in all
!     directions belonging to octant number "oct" is set to outgoing
!     flux using the list of corresponding reflected directions "rdir".
!     "side" takes value of 1 for minimum and 2 for maximum coordinate
!     of boundary plane.
!     "nb", "ns" and "ndir" are also in input and represent respectively
!     number of space moments, number of surfaces on specified side and
!     number of directions per octant.

      REAL     bflx(nb*ns,2,*)
      INTEGER  side,oct,rdir(*)
      INTEGER  d,s

      DO d=(oct-1)*ndir+1,oct*ndir
         DO s=1,nb*ns
            bflx(s,side,d)=bflx(s,side,rdir(d))
         ENDDO
      ENDDO

      END
