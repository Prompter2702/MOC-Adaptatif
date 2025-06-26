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
      REAL, PARAMETER    :: tolcor   = 1.0e-3
    !   REAL, PARAMETER    :: tolcor   = -1.0
      REAL, PARAMETER    :: delt = 10.0
      INTEGER, PARAMETER :: nmat = 2, nsurf = 6
      INTEGER, PARAMETER :: maxinner = 100
      INTEGER, PARAMETER :: n8=8, ns=4
      REAL, PARAMETER    :: rmaj = 2.0, rmin = 1.0 ! Major and minor radii torus

      INTEGER, PARAMETER :: in=1,out=2

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

      INTEGER, ALLOCATABLE :: x_pixel(:,:) ! list of pixel surfaces in X direction
      INTEGER, ALLOCATABLE :: y_pixel(:,:) ! list of pixel surfaces in Y direction
      INTEGER, ALLOCATABLE :: z_pixel(:,:)
      INTEGER, ALLOCATABLE :: iosurf_map(:,:,:)
      INTEGER, ALLOCATABLE :: list_size_reg(:)

      REAL, ALLOCATABLE    :: asrc(:,:,:,:), asrc0(:,:)
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
      INTEGER :: lastadd,cnt,g,i,j,k,nb_cell, nb_pix_tor
      LOGICAL :: ok,oksrc

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      REAL, ALLOCATABLE  :: ccof(:,:,:),icof(:,:,:)
      REAL, ALLOCATABLE  :: ecof(:,:,:),tcof(:,:,:)
      REAL, ALLOCATABLE  :: ccof8(:,:,:,:),icof8(:,:,:,:)
      REAL, ALLOCATABLE  :: ecof8(:,:,:,:),tcof8(:,:,:,:)
      REAL, ALLOCATABLE  :: xstlvl1(:,:)
      REAL(KIND=8),ALLOCATABLE  :: errcor(:,:),errmul(:,:)
      REAL :: errbnd(3)

      REAL(KIND=8), ALLOCATABLE :: Cm(:,:,: ,:,:,:)
      REAL(KIND=8), ALLOCATABLE :: Em(:,:,: ,:,:,:)
      REAL(KIND=8), ALLOCATABLE :: Im(:,:,: ,:,:,:)
      REAL(KIND=8), ALLOCATABLE :: Tm(:,:,: ,:,:,:)


      !Variables à partir des paramètres
      ng = 1
      nc = 4
      nb = 3
      ndim = 3
      nx = 16
      ny = 16
      nz = 16
      nani=0
      nhrm=1 
      nh = 1

      nreg = nmat

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
      ALLOCATE(aflx(nn,nr,nc,noct))

      
      ALLOCATE(ccof(nn,nc,nc),icof(nn,nc,nbd))
      ALLOCATE(ecof(nn,nbd,nc),tcof(nn,nbd,nbd))
      ALLOCATE(ccof8(nn,nc,nc,n8),icof8(nn,nc,nbd,n8))
      ALLOCATE(ecof8(nn,nbd,nc,n8),tcof8(nn,nbd,nbd,n8))
      ALLOCATE(xstlvl1(ng,n8))
      ALLOCATE(errcor(ng, ndir),errmul(ng, ndir))    


      ALLOCATE(sigt(ng, nmat), sigs(ng,0:nani,nmat))
      ALLOCATE(sphr(nhrm,nd))
      ALLOCATE(mu(ndir),eta(ndir),ksi(ndir),w(ndir))
      ALLOCATE(sgnc(nc,nc,noct),sgni(nc,nbd,noct))
      ALLOCATE(sgne(nc,nbd,noct),sgnt(nbd,nbd,noct))
      ALLOCATE(zreg(nr))
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

      ALLOCATE(x_pixel(nbfx,2)) ! list of pixel surfaces in X direction
      ALLOCATE(y_pixel(nbfy,2)) ! list of pixel surfaces in Y direction
      ALLOCATE(z_pixel(nbfz,2)) ! list of pixel surfaces in Z direction
      ALLOCATE(iosurf_map(nsurf,2,noct))
      ALLOCATE(list_size_reg(nmat))


      ALLOCATE(Cm(nn,nc,nreg,nc,nreg,noct)) 
      ALLOCATE(Em(nn,nb,nsurf/2,nc,nreg,noct)) 
      ALLOCATE(Im(nn,nc,nreg,nb,nsurf/2,noct)) 
      ALLOCATE(Tm(nn,nb,nsurf/2,nb,nsurf/2,noct))

      ALLOCATE(asrc0(nn,nc))



      CALL SYSTEM_CLOCK(count_start, count_rate)  ! Capture début
  
      delt3 = (/delt, delt, delt/)
      lastadd = 2

      ! Change nmat according to the number of materials
      ! Definition of the cross sections

    !   CALL sphere_voxel(nx,ny,nz,rmaj,delt/nx,zreg)


      ! géométrie et attribution des régions

      CALL torus_voxel(nx,ny,nz,rmaj,rmin,delt/nx,nb_pix_tor,zreg)
      list_size_reg(1) = nr-nb_pix_tor
      list_size_reg(2) = nb_pix_tor

      ! Création des surfaces 

      x_pixel(:,1) = 1
      x_pixel(:,2) = 2
      y_pixel(:,1) = 3
      y_pixel(:,2) = 4
      z_pixel(:,1) = 5
      z_pixel(:,2) = 6

      ! Création de iosurf

      DO oct=1,noct

        iosurf_map(1,1,oct) = 1
        iosurf_map(2,1,oct) = 1
        
        IF (xinc(oct)==1) THEN
            iosurf_map(1,2,oct) = 1
            iosurf_map(2,2,oct) = 2
        ELSE
            iosurf_map(1,2,oct) = 2
            iosurf_map(2,2,oct) = 1
        ENDIF

        iosurf_map(3,1,oct) = 2
        iosurf_map(4,1,oct) = 2
        IF (yinc(oct)==1) THEN
            iosurf_map(3,2,oct) = 1
            iosurf_map(4,2,oct) = 2
        ELSE
            iosurf_map(3,2,oct) = 2
            iosurf_map(4,2,oct) = 1
        ENDIF

        iosurf_map(5,1,oct) = 3
        iosurf_map(6,1,oct) = 3
        IF (zinc(oct)==1) THEN
            iosurf_map(5,2,oct) = 1
            iosurf_map(6,2,oct) = 2
        ELSE
            iosurf_map(5,2,oct) = 2
            iosurf_map(6,2,oct) = 1
        ENDIF

      ENDDO


      sigt(:,1)   = 0.3
      sigs(:,0,1) = 0.0
      sigt(:,2)   = 0.55
      sigs(:,0,2) = 0.0

      ! Initials Boundary conditions
      bflx = 0.0
      bfly = 0.0
      bflz = 0.0
      
      ! Initialisation of the angular flux and flux moments
      aflx = 0.0 
      flxm = 0.0

      ! Definition of the external source term
      srcm = 0.0

        
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


    
      CALL HCC_COEFF(ng,nn,nmat,nsurf,nsurf/2,nsurf/2,
     &               nc,nb,1,x_pixel,y_pixel,z_pixel,
     &               zreg,iosurf_map,
     &               list_size_reg,
     &               nr,nh,nmat,
     &               nbd,nbfx,nbfy,nbfz,
     &               nx,ny,nz,
     &               nani,nhrm,nd,ndir,noct,
     &               sigt,sigs,
     &               asrc(1,1,1,1),aflx(1,1,1,1),bflx(1,1,1,1,1),
     &               bfly(1,1,1,1,1), bflz(1,1,1,1,1),
     &               mu,eta,ksi,w,pisn,sphr,
     &               zreg, sgnc,sgni, sgne,sgnt,
     &               delt3,pdslu4, aflxmean,
     &               tolcor,asrc0,
     &               aflx0, aflx1,
     &               asrcm0, asrcm1,
     &               finc1, fout1,
     &               finc0, fout0,
     &               ccof,icof,
     &               ecof,tcof,
     &               ccof8,icof8,
     &               ecof8,tcof8, 
     &               Cm,Em,Im,Tm, errcor,errmul)


      print *, "uboundcmext", ubound(Cm)


      print*, "Cm", MAXVAL(ABS(Cm))
      print*, "Em", MAXVAL(ABS(Em))
      print*, "Im", MAXVAL(ABS(Im))
      print*, "Tm", MAXVAL(ABS(Tm))



    !   print *, "Em", Em
    !   print *, "Im", Im
    !   print *, "Tm", Tm


      print *, "HCC_coeff terminé"


      oct = 1
      ok = .FALSE.
      oksrc = .FALSE.


      CALL SYSTEM_CLOCK(count_end, count_rate)    ! Capture fin
      PRINT *, "Temps:", REAL(count_end - count_start)/count_rate,
     &                 " secondes"


      ! Affichage des résultats

      ! Total mean
    !   print*,"Moyenne flux", SUM(flxm(1,:,1,1,lastadd), dim=1)/nr
    !   print*,"Moyenne flux", SUM(aflx(1,:,1,1), dim=1)/nr


    !   flxmean = 0.0
    !   DO oct=1,noct
    !     fst = 0
    !     da = (oct-1)*ndir+1
    !     DO d=1,ndir
    !         flxmean(:ng,:nr) = flxmean(:ng,:nr) +
    !  &          (sphr(1,da+d-1)*w(d))*aflxmean(fst+1:fst+ng,:nr,oct)
    !         fst = fst + ng
    !     ENDDO
    !   ENDDO 

    !   print*,"Moyenne flux", SUM(flxmean(1,:), dim=1)/nr


      ! Enregistrement des résultats dans un fichier VTK


      WRITE(name, '(A,I0,A,I0,A,I0,A)') "milieux.vtk"
      CALL VOLVTK(nx,ny,nz,
     &            .TRUE.,.TRUE.,.TRUE.,
     &            0.0,0.0,0.0,                     
     &            (/delt, delt, delt/),
     &            REAL(zreg),name)

    !   WRITE(name, '(A,I0,A,I0,A,I0,A)') "flx_vol.vtk"
    !   CALL VOLVTK(nx,ny,nz,
    !  &            .TRUE.,.TRUE.,.TRUE.,
    !  &            0.0,0.0,0.0,                     
    !  &            (/delt, delt, delt/),
    !  &            aflx(1,:,1,1),name)

      print *, 'FIN'
    !  
      END PROGRAM