!****************************************************************************
!     "nb"    is nubmer of spatial moments for mesh boundary flux.
!     "nc"    is number of spatial moments for mesh interior flux.
!     "nbd"   is a dimension of the coefficient matrices nb*nd

      PROGRAM SWEEP_ADAPTATIV

      USE SNQHRM
      USE SWEEP8ONE

      IMPLICIT NONE

      INTEGER, PARAMETER :: ng = 1, ndir = 1, ndim = 3
      INTEGER, PARAMETER :: nn = ndir*ng
      INTEGER, PARAMETER :: nc=4, nb=3, nd = ndir*8
      INTEGER, PARAMETER :: nbd = nb*ndim
      INTEGER, PARAMETER :: nx = 4
      INTEGER, PARAMETER :: ny = 4
      INTEGER, PARAMETER :: nz = 4
      INTEGER, PARAMETER :: nbfx = ny*nz,nbfy=nx*nz,nbfz=nx*ny
      INTEGER, PARAMETER :: nani=0,nhrm=1 
      INTEGER, PARAMETER :: nr = nx*ny*nz , nh = 1
      REAL, PARAMETER :: delt = 1.0
      INTEGER, PARAMETER :: nmat = 1
! Nombre de milieux

      INTEGER, PARAMETER :: noct = 8

      REAL    :: flxm(ng,nr,nh,nc)
      REAL    :: bflx(nn,nb,nbfx,2,noct),bfly(nn,nb,nbfy,2,noct),
     &           bflz(nn,nb,nbfz,2,noct)
      REAL    :: sigt(ng, nmat), sigs(ng,0:nani,nmat)
      REAL(KIND=8) :: mu(ndir),eta(ndir),ksi(ndir),w(ndir),pisn
      REAL    :: sphr(nhrm,nd)
      INTEGER :: sgnc(nc,nc,noct),sgni(nc,nbd,noct)
      INTEGER :: sgne(nc,nbd,noct),sgnt(nbd,nbd,noct)
      INTEGER :: zreg(nr)
      REAL    :: aflx(nn,nr,nc,noct),asrc(nn,nr,nc,noct)
      REAL    :: dsrc(nn,nr,nc,noct)
      REAL    :: finc(nn,nbd),fout(nn,nbd)
      REAL    :: flxp(ng,nr,nh,nc),srcm(ng,nr,nh,nc),tmom(ng,nr,nh,nc)
      REAL    :: sigg(ng,nr)
      REAL    :: xaux(nn,nb),yaux(nn,nx,nb),zaux(nn,nx,ny,nb)
      INTEGER :: rdir(nd),dira(nr),dirf(nr)
      LOGICAL :: lgki =.FALSE.
      REAL    :: delt3(3)

      INTEGER :: count_start, count_end, count_rate
      INTEGER :: x,y,z,r
  
      CALL SYSTEM_CLOCK(count_start, count_rate)  ! Capture début
  
      delt3 = (/delt, delt, delt/)

      zreg = 1
      sigs = 0.0
      sigt = 1.0

    !   bflx = 0.0
    !   bfly = 0.0
    !   bflz = 0.0

      bflx(:,1,:,1,:) = 2.0
      bfly(:,1,:,1,:) = 2.0    
      bflz(:,1,:,1,:) = 2.0

      aflx = 1.0 
      asrc = 0.0
      srcm = 0.0

    !   DO z=nz/2,nz/2+1
    !     DO y=ny/2,ny/2+1
    !         DO x=nx/2,nx/2+1
    !             r= x + (y-1)*nx + (z-1)*nx*ny
    !             srcm(:,r,1,1) = 8.0
    !         ENDDO
    !     ENDDO
    !   ENDDO

    !   DO z=1,2
    !     DO y=1,2
    !         DO x=1,2
    !             r= x + (y-1)*nx + (z-1)*nx*ny
    !             srcm(:,r,1,1) = 12.0
    !         ENDDO
    !     ENDDO
    !   ENDDO

    !   print *,"av flux", aflx(:,:,1,1)
      
      CALL COFSGN(sgnc,sgni,sgne,sgnt,nc,nbd,ndim,2)
      CALL SNQDLFT(ndir,nd,nani,ndim, nhrm, mu,eta,ksi,w,sphr)
      
      CALL SWEE3D_ADAPTIVE(nn,ng,nr,nh,nc,nmat,
     &                     nb,nbd,nbfx,nbfy,nbfz,
     &                     nx,ny,nz,
     &                     nani,nhrm,nd,ndir,
     &                     sigt,sigs,
     &                     flxp,srcm,
     &                     bflx,bfly,bflz,
     &                     mu,eta,ksi,w,pisn,sphr,
     &                     sgnc,sgni,sgne,sgnt,
     &                     finc,fout,tmom,sigg,
     &                     aflx,asrc,dsrc,
     &                     xaux,yaux,zaux,
     &                     rdir,
     &                     zreg, 
     &                     dira,dirf,
     &                     lgki,
     &                     flxm, delt3)

      print *,flxm(1,:,1,1)

      CALL SYSTEM_CLOCK(count_end, count_rate)    ! Capture fin
      PRINT *, "Temps:", REAL(count_end - count_start)/count_rate,
     &                 " secondes"


      print *,"bflx", SUM(bflx(1,1,:,2,1), dim=1)/(ny*nz)
      print *,"bfly", SUM(bfly(1,1,:,2,1), dim=1)/(nx*nz)
      print *,"bflz", SUM(bflz(1,1,:,2,1), dim=1)/(nx*ny)

      print *,"Flux final", aflx(:,:,1,1)

      !   CALL VOLVTK(1,ny,nz,
    !  &           .FALSE.,.TRUE.,.TRUE.,
    !  &            delt, 0.0, 0.0,
    !  &            (/delt, delt, delt/),
    !  &            bflx(1,1,:,2,1),"boundx.vtk")
      
    !   CALL VOLVTK(nx,1,ny,
    !  &           .TRUE.,.FALSE.,.TRUE.,
    !  &             0.0, delt, 0.0,
    !  &            (/delt, delt, delt/),
    !  &            bfly(1,1,:,2,1),"boundy.vtk")
      
    !   CALL VOLVTK(nx,ny,1,
    !  &           .TRUE.,.TRUE.,.FALSE.,
    !  &            0.0, 0.0, delt, 
    !  &            (/delt, delt, delt/),
    !  &            bflz(1,1,:,2,1),"boundz.vtk")
     
    !   CALL VOLVTK(nx,1,nz,0.0,delt,0.0, (/delt, delt, delt/),
    !  &           bfly(1,1,:,2,1),"boundy.vtk")
    !   CALL VOLVTK(nx,ny,1,0.0,0.0,delt, (/delt, delt, delt/),
    !  &            bflz(1,1,:,2,1) ,"boundz.vtk")


       CALL VOLVTK(nx,ny,nz,
     &             .TRUE.,.TRUE.,.TRUE.,
     &              0.0, 0.0, 0.0,                     
     &              (/delt, delt, delt/),
     &              flxm(1,:,1,1),"flxmtt.vtk")

      print *, 'FIN'
    !  
      END PROGRAM