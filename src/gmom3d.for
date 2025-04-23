      SUBROUTINE GMOM3D(ng,na,nh,nc,nr,sigs,flxp,srcm,zreg,sigg,tmom)

!     Adds the contribution of selfscattering sources "sigs*flxp"
!     to source moments "srcm". The result is returned in "tmom".
      IMPLICIT NONE
      INTEGER :: ng,na,nr,nh,nc
      REAL    :: flxp(ng,nr,nh,nc),srcm(ng,nr,nh,nc),tmom(ng,nr,nh,nc)
      REAL    :: sigs(ng,0:na,*),sigg(ng,nr)
      INTEGER :: zreg(*)
      INTEGER :: a,c,r,h,k,m


      DO a=0,na
         DO r=1,nr
            m=zreg(r)
            sigg(:,r)= sigs(:,a,m) ! *tlp1
         ENDDO
         h=a*a+1
         DO c=1,nc
         DO k=0,2*a
            tmom(:,:,h+k,c)=srcm(:,:,h+k,c)+(sigg(:,:)*flxp(:,:,h+k,c))
         ENDDO
         ENDDO
      ENDDO

      END SUBROUTINE
      
      

      SUBROUTINE GIRSRC(nhrm,ng,ndir,nr,nh,nc,sphr,tmom,asrc)

!     Given the source moments "tmom" and spherical harmonics
!     coefficients "sphr" computes the directional source "srcm" as
!
!     asrc(r,c) = asrc(r,c) + Sum_{h} sphr(h)*tmom(r,h,c)
      IMPLICIT NONE 
      INTEGER :: nhrm,ng,ndir,nr,nh,nc
      REAL    :: tmom(ng,nr,nh,nc)
      REAL    :: sphr(nhrm,ndir)
      REAL    :: asrc(ng,ndir,nr,nc)
      INTEGER :: c,h,d

      DO c=1,nc
      DO d =1,ndir 
         asrc(:ng,d,:nr,c)=sphr(1,d)*tmom(:ng,:nr,1,c)
         DO h=2,nh
            asrc(:ng,d,:nr,c) = asrc(:ng,d,:nr,c) + 
     &                          sphr(h,d)*tmom(:ng,:nr,h,c)
         ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE