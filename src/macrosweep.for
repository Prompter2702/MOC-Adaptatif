      SUBROUTINE MACRO_SWEEP3D(nn,ng,nr,nh,nc,nmat,
     &                          ncelx,ncely,ncelz,
     &                          nb,nbd,nx,ny,nz, delt3,
     &                          nani,nhrm,nd,ndir,
     &                          sigt,sigs,
     &                          bflx,bfly,bflz,
     &                          mu,eta,ksi,w,sphr,pisn,
     &                          srcm,sigg,
     &                          dsrc,
     &                          rdir,
     &                          zreg, 
     &                          dira,dirf,
     &                          lgki,
     &                          flxm,
     &                          maxinner, tolinner)

      USE FLGINC
      USE SNQHRM
! nr nbre de pixels par cellules = nx*ny*nz

      IMPLICIT NONE

      INTEGER, PARAMETER :: noct = 8, ndim=3

      INTEGER , INTENT(IN) :: nn,ng,nr,nh,nc,nmat,ncelx,ncely,ncelz
      INTEGER , INTENT(IN) :: nb,nbd,nx,ny,nz,nani,nhrm,nd,ndir
      REAL, INTENT(IN)     :: sigt(ng,nmat),sigs(ng,0:nani,nmat)
      REAL, INTENT(IN)     :: delt3(3)
      REAL, INTENT(IN)     :: tolinner
      INTEGER, INTENT(IN)  :: maxinner

      REAL, INTENT(INOUT) ::bflx(nn,nb,ny*nz,2,noct,ncelx,ncely,ncelz)
      REAL, INTENT(INOUT) ::bfly(nn,nb,nx*nz,2,noct,ncelx,ncely,ncelz)
      REAL, INTENT(INOUT) ::bflz(nn,nb,nx*ny,2,noct,ncelx,ncely,ncelz)
      
      REAL(KIND=8), INTENT(IN) :: mu(ndir),eta(ndir),ksi(ndir), 
     &               w(ndir),sphr(nhrm,nd),pisn
      REAL, INTENT(INOUT)  :: srcm(ng,nr,nh,nc,ncelx,ncely,ncelz)
    !   REAL, INTENT(INOUT)  :: tmom(ng,nr,nh,nc)
      
      REAL, INTENT(INOUT)   :: flxm(ng,nr,nh,nc,ncelx,ncely,ncelz)
      INTEGER, INTENT(IN) :: zreg(nr,ncelx,ncely,ncelz)
      REAL, INTENT(INOUT)    :: sigg(ng,nr)
      
      !   REAL    :: aflx(nn,nr,nc,noct)
      REAL    :: dsrc(nn,nr,nc,*)
      REAL    :: asrc(nn,nr,nc,noct)
      INTEGER :: rdir(nd,*),dira(*),dirf(*)
      LOGICAL :: lgki      

      INTEGER :: sgnc(nc,nc,noct),sgni(nc,nbd,noct)
      INTEGER :: sgne(nbd,nc,noct),sgnt(nbd,nbd,noct)
      
      REAL :: errinner
      INTEGER :: cntinner, oc,oct,x,y,z, xout,yout,zout

      
      CALL COFSGN(sgnc,sgni,sgne,sgnt,nc,nbd,ndim,2)
      
      errinner = tolinner+1.0
      cntinner = 0



      DO WHILE(errinner>tolinner .AND. cntinner<20)
        ! errinner = 0.0
        cntinner = cntinner + 1

        DO z=1,ncelz
        DO y=1,ncely
        DO x=1,ncelx

            ! print *,"bflx", bflx(:,:,:,2,7,x,y,z)

          CALL SWEE3D_ADAPTIVE(nn,ng,nr,nh,nc,nmat,
     &                     nb,nbd,ny*nz,nx*nz,nx*ny,
     &                     nx,ny,nz,
     &                     nani,nhrm,nd,ndir,
     &                     sigt,sigs,
     &                     srcm(:,:,:,:,x,y,z),
     &                     bflx(:,:,:,:,:,x,y,z),
     &                     bfly(:,:,:,:,:,x,y,z),
     &                     bflz(:,:,:,:,:,x,y,z),
     &                     mu,eta,ksi,w,pisn,sphr,
     &                     sgnc,sgni,sgne,sgnt,
     &                     sigg,
     &                     dsrc,
     &                     rdir,
     &                     zreg(:,x,y,z), 
     &                     dira,dirf,
     &                     lgki,
     &                     flxm(:,:,:,:,x,y,z), delt3,
     &                     maxinner,tolinner)

            ! print *,"bflx2", bflx(:,:,:,1,7,x,y,z)
            ! print *,"bflx" ,bflx(:,:,:,:,:,x,y,z)
            ! print *,"bfly" ,bfly(x,y,z,:,:,:,:,:)
            ! print *,"bflz" ,bflz(x,y,z,:,:,:,:,:)
            ! print *, "flxm",flxm(:,:,:,:,x,y,z)

      ENDDO
      ENDDO
      ENDDO


      DO oct=1,noct
        oc=olst3D(oct)
        xout=3-xinc(oc)
        yout=3-yinc(oc)
        zout=3-zinc(oc)

        DO x=1,ncelx-1
          bflx(:,:,:,xinc(oc),oc,x+2-xinc(oc),:,:) = 
     &                          bflx(:,:,:,xout,oc,x+xinc(oc)-1,:,:) 
        ENDDO

        DO y=1,ncely-1
          bfly(:,:,:,yinc(oc),oc,:,y+2-yinc(oc),:) = 
     &                          bfly(:,:,:,yout,oc,:,y+yinc(oc)-1,:)
        ENDDO

        DO z=1,ncelz-1
          bflz(:,:,:,zinc(oc),oc,:,:,z+2-zinc(oc)) =
     &                          bflz(:,:,:,zout,oc,:,:,z+zinc(oc)-1)
        ENDDO

      ENDDO 


      !Iterations
      ENDDO

      print *,"cntinner",cntinner




      END SUBROUTINE MACRO_SWEEP3D