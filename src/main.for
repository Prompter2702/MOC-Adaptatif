!****************************************************************************
!     "nb"    is nubmer of spatial moments for mesh boundary flux.
!     "nc"    is number of spatial moments for mesh interior flux.
!     "nbd"   is a dimension of the coefficient matrices nb*nd

      PROGRAM SWEEP_ADAPTATIV
      USE SNQHRM
      USE SWEEP8ONE
      USE FLGINC
      USE SRCCORR
      USE FLXCOM

      IMPLICIT NONE

      INTEGER, PARAMETER :: n = 2 ! Formule triangle
      INTEGER, PARAMETER :: ityp = 4 ! ChebyshevDoubleLegendre
      REAL, PARAMETER    :: tolinner = 1.0e-6
      REAL, PARAMETER    :: tolcor   = 1.0e-2
    !   REAL, PARAMETER    :: tolcor   = -1.0
      REAL, PARAMETER    :: delt = 1.0
      INTEGER, PARAMETER :: nmat = 2
      INTEGER, PARAMETER :: maxinner = 100
      INTEGER, PARAMETER :: n8=8, ns=4
      INTEGER            :: nn
      

      REAL(KIND=8), ALLOCATABLE :: mu(:),eta(:),ksi(:),w(:)

      REAL,    ALLOCATABLE :: flxm(:,:,:,:,:)
      REAL,    ALLOCATABLE :: bflx(:,:,:,:,:),
     &                        bfly(:,:,:,:,:),
     &                        bflz(:,:,:,:,:)
      REAL,    ALLOCATABLE :: sigt(:,:), sigs(:,:,:)
      REAL,    ALLOCATABLE :: sphr(:,:)
      INTEGER, ALLOCATABLE :: sgnc(:,:,:),sgni(:,:,:)
      INTEGER, ALLOCATABLE :: sgne(:,:,:),sgnt(:,:,:)
      INTEGER, ALLOCATABLE :: zreg(:)
      REAL,    ALLOCATABLE :: aflx(:,:,:,:)
      REAL,    ALLOCATABLE :: dsrc(:,:,:,:)
      REAL,    ALLOCATABLE :: srcm(:,:,:,:),tmom(:,:,:,:)
      REAL,    ALLOCATABLE :: sigg(:,:)
      INTEGER, ALLOCATABLE :: rdir(:,:),dira(:),dirf(:)
      REAL,    ALLOCATABLE :: aflxmean(:,:,:), flxmean(:,:)

      REAL, ALLOCATABLE    :: asrc(:,:,:,:)
      REAL(KIND=8), ALLOCATABLE    :: pdslu4(:)
      REAL, ALLOCATABLE    :: aflx0(:,:), aflx1(:,:,:)
      REAL, ALLOCATABLE    :: asrcm0(:,:,:), asrcm1(:,:,:,:)
      REAL, ALLOCATABLE    :: finc1(:,:,:,:), fout1(:,:,:,:)
      REAL, ALLOCATABLE    :: finc0(:,:,:), fout0(:,:,:)

      REAL    :: pisn
      REAL    :: delt3(3)
      INTEGER :: count_start, count_end, count_rate
      INTEGER :: x,y,z,r,oct, d,fst,da
      CHARACTER(LEN=20) :: name
      INTEGER :: lastadd

    
      !Variables à partir des paramètres
      ng = 1
      nc = 4
      nb = 3
      ndim = 3
      nx = 32
      ny = 32
      nz = 32
      nani=0
      nhrm=1 
      nh = 1

      nbfx = ny*nz
      nbfy = nx*nz
      nbfz = nx*ny
      nbf = (/nbfx,nbfy,nbfz/)
      nr = nx*ny*nz
      ndir = (n/2)*((n/2)+1)/2
      nd = ndir*8
      nn = ndir*ng  
      nbd = nb*ndim
      noct = 2**ndim
      lgki =.FALSE.

      ! ALLOCATION DES TABLEAUX
      ALLOCATE(flxm(ng,nr,nh,nc,2))
      ALLOCATE(bflx(nn,nb,nbfx,2,noct))
      ALLOCATE(bfly(nn,nb,nbfy,2,noct))
      ALLOCATE(bflz(nn,nb,nbfz,2,noct))

      ALLOCATE(sigt(ng, nmat), sigs(ng,0:nani,nmat))
      ALLOCATE(sphr(nhrm,nd))
      ALLOCATE(mu(ndir),eta(ndir),ksi(ndir),w(ndir))
      ALLOCATE(sgnc(nc,nc,noct),sgni(nc,nbd,noct))
      ALLOCATE(sgne(nc,nbd,noct),sgnt(nbd,nbd,noct))
      ALLOCATE(zreg(nr))
      
      ALLOCATE(aflx(nn,nr,nc,noct))
      ALLOCATE(dsrc(nn,nr,nc,noct))
      
      ALLOCATE(srcm(ng,nr,nh,nc),tmom(ng,nr,nh,nc))
      ALLOCATE(sigg(ng,nr))
      ALLOCATE(rdir(nd,3))
      ALLOCATE(dira(nr),dirf(nr))
      ALLOCATE(aflxmean(nn,nr,noct), flxmean(ng,nr))
      ALLOCATE(asrc(nn,nr,nc,noct))
      ALLOCATE(pdslu4(ndir))
      ALLOCATE(aflx0(nn,nc), aflx1(nn,nc,n8))
      ALLOCATE(asrcm0(ng,ndir,nc), asrcm1(ng,ndir,nc,n8))
      ALLOCATE(finc1(nn,nb,3,ns), fout1(nn,nb,3,ns))
      ALLOCATE(finc0(nn,nb,3), fout0(nn,nb,3))



      CALL SYSTEM_CLOCK(count_start, count_rate)  ! Capture début
  
      delt3 = (/delt, delt, delt/)
      lastadd = 2

      ! Change nmat according to the number of materials
      ! Definition of the cross sections
      sigt(:,1)   = 1.0
      sigs(:,0,1) = 0.8
      sigt(:,2)   = 20.0
      sigs(:,0,2) = 0.0
    !   sigt(:,3)   = 0.0
    !   sigs(:,0,3) = 0.0
      
      sigg = 0.0

      ! Geometry settings
      zreg = 1

      DO z = 13,16
      DO y = 13,16
      DO x = 13,16
        r = x + (y-1)*nx + (z-1)*nx*ny
        zreg(r) = 2
      ENDDO
      ENDDO
      ENDDO

      ! Initials Boundary conditions
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
    !   DO z=6,6
    !   DO y=6,6
    !   DO x=6,6
    !     r= x + (y-1)*nx + (z-1)*nx*ny
    !     srcm(:,r,1,1) = 15.0
    !   ENDDO
    !   ENDDO
    !   ENDDO

        
      DO oct=1,noct
        da = (oct-1)*ndir+1
        print *,"oct", oct
        print *, xinc(oct), yinc(oct), zinc(oct)
      ENDDO

      ! Creation of the signe matrices
      CALL COFSGN(sgnc,sgni,sgne,sgnt,nc,nbd,ndim,2)
      ! Creation of the angular quadrature
      CALL SNSETC(n,n,ityp,mu,eta,ksi,w)
      w = w/4
      ! Creating the spherical harmonics
      CALL SNQDLFT(ndir,nd,nani,ndim, nhrm, mu,eta,ksi,w,sphr)

      CALL GAUSS_LU4(100,ndir, delt3, mu,eta,ksi, pdslu4)

    !   print *, mu, eta, ksi


      ! Main function that do the total flux calculation
      CALL SWEE_TOT_3D_ADAPTIVE(nn,ng,nr,nh,nc,nmat,
     &                          nb,nbd,nbfx,nbfy,nbfz,
     &                          nx,ny,nz, 
     &                          nani,nhrm,nd,ndir,
     &                          sigt,sigs,srcm,asrc,
     &                          bflx,bfly,bflz,
     &                          mu,eta,ksi,w,pisn,sphr,
     &                          sgnc,sgni,sgne,sgnt,
     &                          sigg,tmom,
     &                          dsrc,aflx,
     &                          rdir,
     &                          zreg, 
     &                          dira,dirf,
     &                          lgki,
     &                          flxm,delt3,
     &                          maxinner, tolinner, tolcor, aflxmean,
     &                          pdslu4,
     &                          aflx0, aflx1,
     &                          asrcm0, asrcm1,
     &                          finc1, fout1,
     &                          finc0, fout0, lastadd)


      CALL SYSTEM_CLOCK(count_end, count_rate)    ! Capture fin
      PRINT *, "Temps:", REAL(count_end - count_start)/count_rate,
     &                 " secondes"


      ! Affichage des résultats

      ! Total mean
      print*,"Moyenne flux", SUM(flxm(1,:,1,1,lastadd), dim=1)/nr



      flxmean = 0.0
      DO oct=1,noct
        fst = 0
        da = (oct-1)*ndir+1
        DO d=1,ndir
            flxmean(:ng,:nr) = flxmean(:ng,:nr) +
     &          (sphr(1,da+d-1)*w(d))*aflxmean(fst+1:fst+ng,:nr,oct)
            fst = fst + ng
        ENDDO
      ENDDO 

      print*,"Moyenne flux", SUM(flxmean(1,:), dim=1)/nr


      ! Enregistrement des résultats dans un fichier VTK

      WRITE(name, '(A,I0,A,I0,A,I0,A)') "flx_vol.vtk"
      CALL VOLVTK(nx,ny,nz,
     &            .TRUE.,.TRUE.,.TRUE.,
     &            0.0,0.0,0.0,                     
     &            (/delt, delt, delt/),
     &            flxm(1,:,1,1,lastadd),name)

      WRITE(name, '(A,I0,A,I0,A,I0,A)') "flx_mean_vol.vtk"
      CALL VOLVTK(nx,ny,nz,
     &            .TRUE.,.TRUE.,.TRUE.,
     &            0.0,0.0,0.0,                     
     &            (/delt, delt, delt/),
     &            flxmean(1,:),name)



      print *, 'FIN'
    !  
      END PROGRAM


          !   typc = 5 ! SpecularReflection
    !   DO d=1,ndir
    !     ! Directions des reflexions sur une face x
    !     rdir((1-1)*ndir + d,1) = (2-1)*ndir + d
    !     rdir((2-1)*ndir + d,1) = (1-1)*ndir + d
    !     rdir((3-1)*ndir + d,1) = (4-1)*ndir + d
    !     rdir((4-1)*ndir + d,1) = (3-1)*ndir + d
    !     rdir((5-1)*ndir + d,1) = (6-1)*ndir + d
    !     rdir((6-1)*ndir + d,1) = (5-1)*ndir + d
    !     rdir((7-1)*ndir + d,1) = (8-1)*ndir + d
    !     rdir((8-1)*ndir + d,1) = (7-1)*ndir + d
    !     ! Directions des reflexions sur une face y
    !     rdir((1-1)*ndir + d,2) = (4-1)*ndir + d
    !     rdir((4-1)*ndir + d,2) = (1-1)*ndir + d
    !     rdir((2-1)*ndir + d,2) = (3-1)*ndir + d
    !     rdir((3-1)*ndir + d,2) = (2-1)*ndir + d
    !     rdir((5-1)*ndir + d,2) = (8-1)*ndir + d
    !     rdir((8-1)*ndir + d,2) = (5-1)*ndir + d
    !     rdir((6-1)*ndir + d,2) = (7-1)*ndir + d
    !     rdir((7-1)*ndir + d,2) = (6-1)*ndir + d
    !     ! Directions des reflexions sur une face z
    !     rdir((1-1)*ndir + d,3) = (5-1)*ndir + d
    !     rdir((5-1)*ndir + d,3) = (1-1)*ndir + d
    !     rdir((2-1)*ndir + d,3) = (6-1)*ndir + d
    !     rdir((6-1)*ndir + d,3) = (2-1)*ndir + d
    !     rdir((3-1)*ndir + d,3) = (7-1)*ndir + d
    !     rdir((7-1)*ndir + d,3) = (3-1)*ndir + d
    !     rdir((4-1)*ndir + d,3) = (8-1)*ndir + d
    !     rdir((8-1)*ndir + d,3) = (4-1)*ndir + d
    !   ENDDO