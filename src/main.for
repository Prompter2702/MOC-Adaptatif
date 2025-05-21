!****************************************************************************
!     "nb"    is nubmer of spatial moments for mesh boundary flux.
!     "nc"    is number of spatial moments for mesh interior flux.
!     "nbd"   is a dimension of the coefficient matrices nb*nd

      PROGRAM SWEEP_ADAPTATIV

      USE SNQHRM
      USE SWEEP8ONE
      USE FLGINC

      IMPLICIT NONE

      INTEGER, PARAMETER :: n = 2   ! Formule triangle
      INTEGER, PARAMETER :: ityp = 4 ! ChebyshevDoubleLegendre
      INTEGER, PARAMETER :: ng = 1, ndir =(n/2)*((n/2)+1)/2, ndim = 3
      INTEGER, PARAMETER :: nn = ndir*ng
      INTEGER, PARAMETER :: nc=4, nb=3, nd = ndir*8
      INTEGER, PARAMETER :: nbd = nb*ndim
      INTEGER, PARAMETER :: nx = 4
      INTEGER, PARAMETER :: ny = 4
      INTEGER, PARAMETER :: nz = 4
      INTEGER, PARAMETER :: nbfx = ny*nz,nbfy=nx*nz,nbfz=nx*ny
      INTEGER, PARAMETER :: nani=0,nhrm=1 
      INTEGER, PARAMETER :: nr = nx*ny*nz , nh = 1
      REAL, PARAMETER    :: delt = 1.0
      INTEGER, PARAMETER :: nmat = 2
      INTEGER, PARAMETER :: maxinner = 500
      REAL, PARAMETER    :: tolinner = 1.0e-4
      REAL, PARAMETER    :: tolcor   = 1.0e-4
    !   REAL, PARAMETER    :: tolcor   = -1.0
! Nombre de milieux

      INTEGER, PARAMETER :: noct = 8

      REAL    :: flxm(ng,nr,nh,nc)
      REAL    :: bflx(nn,nb,nbfx,2,noct ),
     &           bfly(nn,nb,nbfy,2,noct ),
     &           bflz(nn,nb,nbfz,2,noct )
      REAL    :: sigt(ng, nmat), sigs(ng,0:nani,nmat)
      REAL(KIND=8) :: mu(ndir),eta(ndir),ksi(ndir),w(ndir),pisn
      REAL    :: sphr(nhrm,nd)
      INTEGER :: sgnc(nc,nc,noct),sgni(nc,nbd,noct)
      INTEGER :: sgne(nc,nbd,noct),sgnt(nbd,nbd,noct)
      INTEGER :: zreg(nr)
      REAL    :: aflx(nn,nr,nc,noct)
      REAL    :: dsrc(nn,nr,nc,noct)
      REAL    :: finc(nn,nbd),fout(nn,nbd)
      REAL    :: flxp(ng,nr,nh,nc),srcm(ng,nr,nh,nc ),
     &          tmom(ng,nr,nh,nc)
      REAL    :: sigg(ng,nr)
      INTEGER :: rdir(nd),dira(nr),dirf(nr)
      LOGICAL :: lgki =.FALSE.
      REAL    :: delt3(3)

      INTEGER :: count_start, count_end, count_rate
      INTEGER :: x,y,z,r,oct,cnt
      CHARACTER(LEN=20) :: name


      CALL SYSTEM_CLOCK(count_start, count_rate)  ! Capture début
  
      delt3 = (/delt, delt, delt/)

      sigt(:,1) = 1.0
      sigs(:,0,1) = 0.9999

      sigt(:,2) = 10.0
      sigs(:,0,2) = 9.999
      
      zreg = 1
      DO z = 2,3
        DO y = 2,3
            DO x = 2,3
            r = x + (y-1)*nx + (z-1)*nx*ny
            zreg(r) = 2
          ENDDO
        ENDDO
      ENDDO

      sigg = 0.0

      bflx = 0.0
      bfly = 0.0
      bflz = 0.0
        bflx(:,1,:,:,:) = 10.0
        bfly(:,1,:,:,:) = 10.0
        bflz(:,1,:,:,:) = 10.0

      aflx = 0.0 
      srcm = 0.0
      flxm = 0.0

    !   srcm(:,:,1,1) = 1.0
    !   DO z=2,3
    !     DO y=2,3
    !         DO x=2,3
    !             r= x + (y-1)*nx + (z-1)*nx*ny
    !             srcm(:,r,1,1) = 5.0
    !         ENDDO
    !     ENDDO
    !   ENDDO
      
      CALL COFSGN(sgnc,sgni,sgne,sgnt,nc,nbd,ndim,2)
      CALL SNSETC(n,n,ityp,mu,eta,ksi,w)
      w = w/4
      CALL SNQDLFT(ndir,nd,nani,ndim, nhrm, mu,eta,ksi,w,sphr)

      CALL SWEE3D_ADAPTIVE(nn,ng,nr,nh,nc,nmat,
     &                          nb,nbd,nbfx,nbfy,nbfz,
     &                          nx,ny,nz, 
     &                          nani,nhrm,nd,ndir,
     &                          sigt,sigs,srcm,
     &                          bflx,bfly,bflz,
     &                          mu,eta,ksi,w,pisn,sphr,
     &                          sgnc,sgni,sgne,sgnt,
     &                          sigg,
     &                          dsrc,aflx,
     &                          rdir,
     &                          zreg, 
     &                          dira,dirf,
     &                          lgki,
     &                          flxm,delt3,
     &                          maxinner, tolinner, tolcor)

    !   print *,flxm(1,:,1,1)

      CALL SYSTEM_CLOCK(count_end, count_rate)    ! Capture fin
      PRINT *, "Temps:", REAL(count_end - count_start)/count_rate,
     &                 " secondes"

    
      print*, SUM(flxm(1,:,1,1), dim=1)/nr

      DO x=1,ndir
        print *,"bflx", SUM(bflx(x,1,:,2,1), dim=1)/nbfx
        print *,"bfly", SUM(bfly(x,1,:,2,1), dim=1)/nbfy
        print *,"bflz", SUM(bflz(x,1,:,2,1), dim=1)/nbfz
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
     

      WRITE(name, '(A,I0,A,I0,A,I0,A)') "flxmtt_volum.vtk"


        CALL VOLVTK(nx,ny,nz,
     &             .TRUE.,.TRUE.,.TRUE.,
     &              delt*(x-1), delt*(y-1), delt*(z-1),                     
     &              (/delt, delt, delt/),
     &              flxm(1,:,1,1),name)

        ! print*, flxm(1,:,1,1,x,y,z)


      print *, 'FIN'
    !  
      END PROGRAM