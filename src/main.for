!****************************************************************************
!     "nb"    is nubmer of spatial moments for mesh boundary flux.
!     "nc"    is number of spatial moments for mesh interior flux.
!     "nbd"   is a dimension of the coefficient matrices nb*nd

      PROGRAM SWEEP_ADAPTATIV

      USE SNQHRM
      USE SWEEP8ONE
      USE FLGINC

      IMPLICIT NONE

      INTEGER, PARAMETER :: ng = 1, ndir = 3, ndim = 3
      INTEGER, PARAMETER :: nn = ndir*ng
      INTEGER, PARAMETER :: nc=4, nb=3, nd = ndir*8
      INTEGER, PARAMETER :: nbd = nb*ndim
      INTEGER, PARAMETER :: nx = 2
      INTEGER, PARAMETER :: ny = 2
      INTEGER, PARAMETER :: nz = 2
      INTEGER, PARAMETER :: ncelx = 2, ncely = 2, ncelz = 2
      INTEGER, PARAMETER :: nbfx = ny*nz,nbfy=nx*nz,nbfz=nx*ny
      INTEGER, PARAMETER :: nani=0,nhrm=1 
      INTEGER, PARAMETER :: nr = nx*ny*nz , nh = 1
      REAL, PARAMETER    :: delt = 1.0
      INTEGER, PARAMETER :: nmat = 1
      INTEGER, PARAMETER :: maxinner = 100
      REAL, PARAMETER    :: tolinner = 1.0e-4

! Nombre de milieux

      INTEGER, PARAMETER :: noct = 8

      REAL    :: flxm(ng,nr,nh,nc,ncelx,ncely,ncelz)
      REAL    :: bflx(nn,nb,nbfx,2,noct,ncelx,ncely,ncelz),
     &           bfly(nn,nb,nbfy,2,noct,ncelx,ncely,ncelz),
     &           bflz(nn,nb,nbfz,2,noct,ncelx,ncely,ncelz)
      REAL    :: sigt(ng, nmat), sigs(ng,0:nani,nmat)
      REAL(KIND=8) :: mu(ndir),eta(ndir),ksi(ndir),w(ndir),pisn
      REAL    :: sphr(nhrm,nd)
      INTEGER :: sgnc(nc,nc,noct),sgni(nc,nbd,noct)
      INTEGER :: sgne(nc,nbd,noct),sgnt(nbd,nbd,noct)
      INTEGER :: zreg(ncelx,ncely,ncelz,nr)
      REAL    :: aflx(nn,nr,nc,noct)
      REAL    :: dsrc(nn,nr,nc,noct)
      REAL    :: finc(nn,nbd),fout(nn,nbd)
      REAL    :: flxp(ng,nr,nh,nc),srcm(ng,nr,nh,nc,ncelx,ncely,ncelz),
     &          tmom(ng,nr,nh,nc)
      REAL    :: sigg(ng,nr)
      INTEGER :: rdir(nd),dira(nr),dirf(nr)
      LOGICAL :: lgki =.FALSE.
      REAL    :: delt3(3)

      INTEGER :: count_start, count_end, count_rate
      INTEGER :: x,y,z,r,oct
      CHARACTER(LEN=20) :: name


      CALL SYSTEM_CLOCK(count_start, count_rate)  ! Capture début
  
      delt3 = (/delt, delt, delt/)

      zreg = 1
      sigs = 0.0
      sigt = 1.0
      sigs(:,0,:) = 0.2
      sigg = 0.0

      bflx = 0.0
      bfly = 0.0
      bflz = 0.0
      flxm = 0.0

      DO oct=1,8
        bflx(:,1,:,xinc(oct),oct, :,:,:) = 0.0
        bfly(:,1,:,yinc(oct),oct, :,:,:) = 0.0    
        bflz(:,1,:,zinc(oct),oct, :,:,:) = 0.0
      END DO

      aflx = 1.0 
      srcm = 0.0

    !   DO z=nz/2,nz/2+1
    !     DO y=ny/2,ny/2+1
    !         DO x=nx/2,nx/2+1
    !             r= x + (y-1)*nx + (z-1)*nx*ny
    !             srcm(:,r,1,1) = 8.0
    !         ENDDO
    !     ENDDO
    !   ENDDO

      DO z=1,nz
      DO y=1,ny
      DO x=1,nx
            r= x + (y-1)*nx + (z-1)*nx*ny
            srcm(:,r,1,1, 1,1,1) = 12.0
          ENDDO
        ENDDO
      ENDDO

    !   print *,"av flux", aflx(:,:,1,1)
      
    !   CALL COFSGN(sgnc,sgni,sgne,sgnt,nc,nbd,ndim,2)
      CALL SNQDLFT(ndir,nd,nani,ndim, nhrm, mu,eta,ksi,w,sphr)

      CALL MACRO_SWEEP3D(nn,ng,nr,nh,nc,nmat,
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

    !   print *,flxm(1,:,1,1)

      CALL SYSTEM_CLOCK(count_end, count_rate)    ! Capture fin
      PRINT *, "Temps:", REAL(count_end - count_start)/count_rate,
     &                 " secondes"


      DO x=1,ndir
        print *,"bflx", SUM(bflx(x,1,:,2,1, 1,1,1), dim=1)/(ny*nz)
        print *,"bfly", SUM(bfly(x,1,:,2,1, 1,1,1), dim=1)/(nx*nz)
        print *,"bflz", SUM(bflz(x,1,:,2,1, 1,1,1), dim=1)/(nx*ny)
      ENDDO


    !   print *,"bflx", bflx(1,1,:,2,1, 1,1,1)
    !   print *,"bflx", bflx(1,1,:,1,1, 2,1,1)

    !   print *,"Flux final", aflx(:,:,1,1)

!!!!! Création des fichiers VTK pour lire dans paraview

    !     CALL VOLVTK(1,ny,nz,
    !  &           .FALSE.,.TRUE.,.TRUE.,
    !  &            delt, 0.0, 0.0,
    !  &            (/delt, delt, delt/),
    !  &            bflx(1,1,:,2,1,),"boundx.vtk")
      
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
     

      DO z=1,ncelz
      DO y=1,ncely
      DO x=1,ncelx
        WRITE(name, '(A,I0,A,I0,A,I0,A)') "flxmtt_",x,"_",y,"_",z,".vtk"


        CALL VOLVTK(nx,ny,nz,
     &             .TRUE.,.TRUE.,.TRUE.,
     &              delt*(x-1), delt*(y-1), delt*(z-1),                     
     &              (/delt, delt, delt/),
     &              flxm(1,:,1,1,x,y,z),name)

        print*, flxm(1,:,1,1,x,y,z)

      ENDDO
      ENDDO
      ENDDO

      print *, 'FIN'
    !  
      END PROGRAM