!****************************************************************************
!     "nb"    is nubmer of spatial moments for mesh boundary flux.
!     "nc"    is number of spatial moments for mesh interior flux.
!     "nbd"   is a dimension of the coefficient matrices nb*nd

      PROGRAM SWEEP_ADAPTATIV

      USE SNQHRM
      USE SWEEP8ONE
      USE FLGINC
      USE SRCCORR

      IMPLICIT NONE

      INTEGER, PARAMETER :: n = 2  ! Formule triangle
      INTEGER, PARAMETER :: ityp = 4 ! ChebyshevDoubleLegendre
      INTEGER, PARAMETER :: ng = 1, ndir =(n/2)*((n/2)+1)/2, ndim = 3
      INTEGER, PARAMETER :: nn = ndir*ng
      INTEGER, PARAMETER :: nc=4, nb=3, nd = ndir*8
      INTEGER, PARAMETER :: nbd = nb*ndim
      INTEGER, PARAMETER :: nx = 16
      INTEGER, PARAMETER :: ny = 16
      INTEGER, PARAMETER :: nz = 16
      INTEGER, PARAMETER :: nbfx = ny*nz,nbfy=nx*nz,nbfz=nx*ny
      INTEGER, PARAMETER :: nani=0,nhrm=1 
      INTEGER, PARAMETER :: nr = nx*ny*nz , nh = 1
      REAL, PARAMETER    :: delt = 1.0
      INTEGER, PARAMETER :: nmat = 2
      INTEGER, PARAMETER :: maxinner = 80
      REAL, PARAMETER    :: tolinner = 1.0e-4
      REAL, PARAMETER    :: tolcor   = 1.0e-4
    !   REAL, PARAMETER    :: tolcor   = -1.0
! Nombre de milieux

      INTEGER, PARAMETER :: noct = 8

      REAL    :: flxm(ng,nr,nh,nc)
      REAL    :: bflx(nn,nb,nbfx,2,noct),
     &           bfly(nn,nb,nbfy,2,noct),
     &           bflz(nn,nb,nbfz,2,noct)
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

      REAL :: fout0(nn,nb,3)
      REAL :: errbnd
      REAL :: tcof(ng,nbd,nbd)


      CALL SYSTEM_CLOCK(count_start, count_rate)  ! Capture début
  
      delt3 = (/delt, delt, delt/)

      ! Change nmat according to the number of materials
      ! Definition of the cross sections
      sigt(:,1)   = 1.0
      sigs(:,0,1) = 0.3
      sigt(:,2)   = 20.0
      sigs(:,0,2) = 0.0
      
      sigg = 0.0

      ! Geometry settings
      zreg = 1
      DO z = 7,8
      DO y = 7,8
      DO x = 7,8
            r = x + (y-1)*nx + (z-1)*nx*ny
            zreg(r) = 2
      ENDDO
      ENDDO
      ENDDO

      ! Boundary conditions
      bflx = 0.0
      bfly = 0.0
      bflz = 0.0
      bflx(:,1,:,1,1) = 12.0
      bfly(:,1,:,1,1) = 12.0
      bflz(:,1,:,1,1) = 12.0

      ! Initialisation of the angular flux and flux moments
      aflx = 0.0 
      flxm = 0.0

      ! Definition of the external source term
      srcm = 0.0
    !   srcm(:,:,1,1) = 1.0
    !   DO z=2,3
    !     DO y=2,3
    !         DO x=2,3
    !             r= x + (y-1)*nx + (z-1)*nx*ny
    !             srcm(:,r,1,1) = 5.0
    !         ENDDO
    !     ENDDO
    !   ENDDO
      
      ! Creation of the signe matrices
      CALL COFSGN(sgnc,sgni,sgne,sgnt,nc,nbd,ndim,2)
      ! Creation of the angular quadrature
      CALL SNSETC(n,n,ityp,mu,eta,ksi,w)
      w = w/4
      ! Creating the spherical harmonics
      CALL SNQDLFT(ndir,nd,nani,ndim, nhrm, mu,eta,ksi,w,sphr)


      ! Main function that do the total flux calculation
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



      CALL SYSTEM_CLOCK(count_end, count_rate)    ! Capture fin
      PRINT *, "Temps:", REAL(count_end - count_start)/count_rate,
     &                 " secondes"


      ! Affichage des résultats
    
      ! Total mean
      print*, SUM(flxm(1,:,1,1), dim=1)/nr

      DO x=1,ndir
        print *,"bflx", SUM(bflx(x,1,:,2,1), dim=1)/nbfx
        print *,"bfly", SUM(bfly(x,1,:,2,1), dim=1)/nbfy
        print *,"bflz", SUM(bflz(x,1,:,2,1), dim=1)/nbfz
      ENDDO


      ! Enregistrement des résultats dans un fichier VTK

      WRITE(name, '(A,I0,A,I0,A,I0,A)') "flx_mean_volu.vtk"
      CALL VOLVTK(nx,ny,nz,
     &            .TRUE.,.TRUE.,.TRUE.,
     &            0.0,0.0,0.0,                     
     &            (/delt, delt, delt/),
     &            flxm(1,:,1,1),name)


      print *, 'FIN'
    !  
      END PROGRAM