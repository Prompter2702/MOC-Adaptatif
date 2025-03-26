      SUBROUTINE TRANP4(bflx,bfly,nb,ns,ndir,side,oct,rdir)

!     Treats Pi/2 rotation condition of a domain having Pi/2 rotational
!     symmetry on one of the axes and translation condition on the same
!     domain on one of the outer boundaries.
!     Incoming angular flux "bflx" for a specified side "side" in all
!     directions belonging to octant number "oct" is set to outgoing
!     flux "bfly" using the list of corresponding exiting directions
!     "rdir".
!     "side" takes value of 1 for minimum and 2 for maximum coordinate
!     of boundary plane.
!     "nb", "ns" and "ndir" are also in input and represent respectively
!     number of space moments, number of surfaces on specified side and
!     number of directions per octant.

      REAL     bflx(nb*ns,2,*),bfly(nb*ns,2,*)
      INTEGER  side,oct,rdir(*)
      INTEGER  d,s

      DO d=(oct-1)*ndir+1,oct*ndir
                ! write (9,*)' side,oct,d,rdir ', side,oct,d,'<-',rdir(d)
         DO s=1,nb*ns
            bflx(s,side,d)=bfly(s,side,rdir(d))
         ENDDO
      ENDDO

      END
      
!------------------------------------------------------------------------      
!------------------------------------------------------------------------      
      
      SUBROUTINE TRANP4_TST(bflx,bfly,nb,nsx,nsy,ndir,side,oct,rdir,
     &                      adx,surx,matx)

!     Treats Pi/2 rotation condition of a domain having Pi/2 rotational
!     symmetry on one of the axes and translation condition on the same
!     domain on one of the outer boundaries. 
!     The mesh on incoming and outcgoing axis is not necessary the same. 
!     Incoming angular fl'ux "bflx" for a specified side "side" in all
!     directions belonging to octant number "oct" is set to outgoing
!     flux "bfly" using the list of corresponding exiting directions
!     "rdir".
!     "side" takes value of 1 for minimum and 2 for maximum coordinate
!     of boundary plane.
!     "nb", "ns" and "ndir" are also in input and represent respectively
!     number of space moments, number of surfaces on specified side and
!     number of directions per octant.
!     "adx","surx", and "matx" are respectively the address, 
!     the surface-to-surface mapping the iterpolaion matrix for 
!     non-conforming incoming/outgoing mashes. 

      REAL     bflx(nb,nsx,2,*),bfly(nb,nsy,2,*)
      INTEGER  side,oct,rdir(*)
      INTEGER  adx(*),surx(*)
      REAL(KIND=8) matx(nb,nb,*)
      ! locals 
      REAL(KIND=8) fla(nb)
      INTEGER  d,s,d1,d2,s1,s2,ss 
      
      d1= (oct-1)*ndir+1 ; d2 = d1+ndir-1
      DO d=d1,d2
         DO s=1,nsx
            fla(:nb) = 0.0D0  ; s1 = adx(s) ; s2 = adx(s+1)-1
            DO ss = s1,s2 
                fla(:nb) =  fla(:nb) + 
     &          MATMUL(matx(:nb,:nb,ss),bfly(:nb,surx(ss),side,rdir(d)))
            ENDDO 
            bflx(:nb,s,side,d)=fla(:nb)
         ENDDO
      ENDDO
      
      END
