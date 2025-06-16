       SUBROUTINE HCC_COFHCC_PIXEL(
    ! inputs HCC main dimensions 
     &                       ng,nr,nsur,nsui,nsuo,
     &                       nc,nb,nphy,
     &                       sur_centroid,
     &                       reg_centroid,
     &                       pixel_to_cmpregion,
     &                       list_size_reg,
     &                       list_pix,
        ! inputs: dimensions for the pixel_grif
     &                       npix,nh,nmat,
     &                       nbd,nbfx,nbfy,nbfz,
     &                       nx,ny,nz,
     &                       nani,nhrm,nd,ndir,noct,
       ! inputs: cross section  
     &                       sigt,sigs,
       ! input/output: boundary flux 
     &                       bflx,bfly,bflz,
       ! input: angular quadrature & spherical harmonics
     &                       mu,eta,ksi,w,pisn,sphr,
       ! input: sing matrices to adapt coefficients to octants 
     &                       sgnc,sgni,sgne,sgnt,
       ! input: auxiliar memory for boundary angular fluxes, angular 
     &                        asrc,srcm,
       ! input: pixel-to-medium array 
     &                       zreg, 
       ! output: angular flux moments 
     &                       aflx,
       ! Delta on each direction of the region
     &                       delt3,
      ! max inner iterations and tolerance for inner iterations
     &                       tolcor,
      ! outputs  
     &                       Cm,Im,Em,Tm )

      USE FLGINC
      
      IMPLICIT NONE

      INTEGER, PARAMETER :: n8 = 8, ns4 = 4, ndim = 3
      
      ! HCC data inputs 
      INTEGER :: ng      ! n. of groups 
      INTEGER :: nr      ! n. of regions in HCC
      INTEGER :: nsur    ! n. of surfaces of the HCC
      INTEGER :: nsui    ! n. of incoming surfaces of the HCC
      INTEGER :: nsuo    ! n. of outgoing surfaces of the HCC
      INTEGER :: nc      ! n. of spatial moments per HCC's volume
      INTEGER :: nb      ! n. of spatial moments per HCC's surface 
      INTEGER :: nphy    ! n. of physical HCCs sharing the same HCC geometry

      ! data for pixels in the HCC (local data)
      
      INTEGER :: npix    ! n. of pixels -> nr 
      INTEGER :: nbd     ! nb * (ndim = 3) 
      INTEGER :: nmat    ! n. of media 
      INTEGER :: nbfx    ! n. of pixel surfaces in X direction
      INTEGER :: nbfy    ! n. of pixel surfaces in Y direction
      INTEGER :: nbfz    ! n. of pixel surfaces in Z direction
      INTEGER :: nx      ! n. of pixels in X direction
      INTEGER :: ny      ! n. of pixels in Y direction
      INTEGER :: nz      ! n. of pixels in Z direction

      INTEGER :: list_size_reg(nr)  ! list of nb of pixel in each region
      INTEGER :: list_pix(npix,3,nr)  ! list of pixels in each region 
      INTEGER :: list_cube_reg(2,3,nr) ! list of coordinates of the circonsrit 
      ! cube of each region
      INTEGER :: x_pixel(nbfx,2) ! list of pixel surfaces in X direction
      INTEGER :: y_pixel(nbfy,2) ! list of pixel surfaces in Y direction
      INTEGER :: z_pixel(nbfz,2) ! list of pixel surfaces in Z direction

      INTEGER :: list_oct_surf(8,nsur) ! list of surfaces per incomming octant

      INTEGER :: oct      ! octant index
      REAL    , DIMENSION(ng,nmat)    :: sigt   ! total XS per HCC's region
      REAL    , DIMENSION(ng,nphy,nr) :: sigs   ! total XS per HCC's region
     
      INTEGER , DIMENSION(nsur) :: iosurf_map  ! incoming/outgoing surface index per HCC's surface (per octant )
       
      REAL(KIND=8) , DIMENSION(2,nsur) :: sur_centroid  
      REAL(KIND=8) , DIMENSION(3,nr)   :: reg_centroid
      REAL, INTENT(IN)                 :: delt3(3)  ! sides' lengths of the HCC's box 
      
      INTEGER :: pixel_to_cmpregion(npix) ! mapping pixel-to-computational region 
      INTEGER :: list_square_surf(2,2,nsur) ! square circonsrcit of each surface 

        
      ! outputs 
      REAL(KIND=8),INTENT(INOUT), DIMENSION(ng,ndir,nc,nr,nc,nr)   :: Cm
      REAL(KIND=8),INTENT(INOUT), DIMENSION(ng,ndir,nb,nsuo,nc,nr) :: Em
      REAL(KIND=8),INTENT(INOUT), DIMENSION(ng,ndir,nc,nr,nb,nsui) :: Im
      REAL(KIND=8),INTENT(INOUT), DIMENSION(ng,ndir,nb,nsuo,nb,nsui)::Tm
      
      !! Que des tableaux auxiliaires pour construire les matrices

      REAL, ALLOCATABLE :: asrc(:,:,:,:) ! spatial moments of the angular source 
      REAL, ALLOCATABLE :: aflx(:,:,:,:) ! angular flux along the Z-normal surface 
      REAL, ALLOCATABLE :: bflx(:,:,:,:)   ! boundary angular flux along the X-normal surface 
      REAL, ALLOCATABLE :: bfly(:,:,:,:)   ! boundary angular flux along the Y-normal surface 
      REAL, ALLOCATABLE :: bflz(:,:,:,:)   ! boundary angular flux along the Z-normal surface 
      REAL, ALLOCATABLE :: asrc0(:,:,:)

      REAL, INTENT(INOUT) :: tmom(ng,ndir,nc,npix), ! source moments
     &                       srcm(ng,npix,nh,nc) ! 
      INTEGER, INTENT(IN) :: zreg(nphy,npix) ! pixel-to-medium assignment      
      
      ! pixel cross sections bufer 
      REAL :: sigg(ng,npix,nphy)
      
      ! angular data
      INTEGER :: noct    ! n. of octants 
      INTEGER :: nh      ! n. of angular moments of the flux 
      INTEGER :: nani    ! max anisotropy order 
      INTEGER :: nhrm    ! n. of spherical harmonics 
      INTEGER :: nd      ! n. of directions
      INTEGER :: ndir    ! n. of directions per octant
      
      REAL(KIND=8), INTENT(IN) :: mu(ndir),eta(ndir),ksi(ndir),
     &                            w(ndir),pisn
      REAL                     :: sphr(nhrm,nd)
      
      !  sign matrices 
      INTEGER,INTENT(IN) :: sgnc(nc,nc,noct),sgni(nc,nbd,noct)
      INTEGER,INTENT(IN) :: sgne(nbd,nc,noct),sgnt(nbd,nbd,noct)
      REAL, INTENT(IN) :: tolcor
      


      REAL  :: asrcreg(ng,ndir,nr,nc)
      
      ! locals 
      INTEGER :: npvol, xin,yin,zin, xout,yout,zout, rin,c,b, rout,i,
     &           sout, cout, bout, sin, lastadd, nn

      LOGICAL :: lgki

      REAL(KIND=8), ALLOCATABLE    :: pdslu4(:)
      nn = ng*ndir
      lgki = .FALSE.
      tolcor = 1.0e-3

      ALLOCATE(pdslu4(ndir))
      ALLOCATE(asrc(ng,ndir,npix,nc)) ! spatial moments of the angular source 
      ALLOCATE(aflx(ng,ndir,npix,nc)) ! angular flux along the Z-normal surface 
      ALLOCATE(bflx(ng,ndir,nb,nbfx))   ! boundary angular flux along the X-normal surface 
      ALLOCATE(bfly(ng,ndir,nb,nbfy))   ! boundary angular flux along the Y-normal surface 
      ALLOCATE(bflz(ng,ndir,nb,nbfz))   ! boundary angular flux along the Z-norma
      ALLOCATE(asrc0(ng,ndir,nc))

      
      xout = 3-xin
      yin = yinc(oct)
      yout = 3-yin 
      zin = zinc(oct)
      zout = 3-zin 

      Cm = 0
      Em = 0
      Im = 0
      Tm = 0

      CALL GAUSS_LU4(100,ndir, delt3, mu,eta,ksi, pdslu4)

      ! n. of volume-source problems 
      
      DO rin=1,nr
      DO c=1,nc 
        ! 1- contruct the volume-source for the couple of indexes (r,c) 
        ! r = HCC region, c = spatial volume component 
        ! 2- projection: project the volume HCC source onto the pixel support  
        asrc = 0
        bflx = 0
        bfly = 0
        bflz = 0
        tmom = 0
        srcm = 0
        asrc0 = 0
        asrc0(:,:,c) = 1.0

        CALL SPLITVOLHCC(ng*ndir,npix,nx,ny,
     &                  asrc0,asrc,reg_centroid(:,rin),
     &                  list_pix(:,:,rin))

    ! 3- solve the pixel volume source problem 
    ! Construire une fonction plus adaptée !!

       CALL SWEEP_ADA_ONE_OCTANT(nn,ng,ndir,nr,nh,
     &                           nx,ny,nz, delt3,
     &                           xinc,yinc,zinc, oct,
     &                           sigt,mu,eta,ksi,w,pdslu4,
     &                           sgnc,sgni,sgne,sgnt,
     &                           zreg,asrc,aflx,
     &                           bflx, bfly, bflz,tolcor)
     
        ! 4- fill-in the C and E matrix column
        CALL MERGEVOLHCC(nn,npix,nr,nx,ny,nz,
     &                   list_cube_reg,pixel_to_cmpregion,
     &                   aflx, Cm(1,1,1,1,c,rin))
        
        CALL MERGEBOUNDHCC(nn,nsuo,nx,ny,nz,
     &                      nbfx, nbfy, nbfz, list_square_surf,
     &                      x_pixel(:,xout), 
     &                      y_pixel(:,yout),
     &                      z_pixel(:,zout),
     &                      list_oct_surf,   
     &                      bflx(1,1,1,1),
     &                      bfly(1,1,1,1),
     &                      bflz(1,1,1,1),
     &                      Em(1,1,1,1,c,rin))
        END DO

      ENDDO 
          
      ! n. of surface-source problems 

      DO sin=1,nsur
      DO b=1,nb 
         ! 1- contruct boundary-source for the couple of indexes (s,b) 
         ! s = HCC surface, b = spatial surface component 
      
         ! 2- projection: project the boundary HCC source onto
         ! the pixel boundary support  
         asrc = 0
         bflx = 0
         bfly = 0
         bflz = 0
         tmom = 0
         srcm = 0
         flxm = 0  
        
         DO i=1,nbfx 
            IF( x_pixel(i,1) == sin) bflx(:,:,b,i,1,:) = 1
            IF( x_pixel(i,2) == sin) bflx(:,:,b,i,2,:) = 1 
         ENDDO
         DO i=1,nbfy 
            IF( y_pixel(i,1) == sin) bfly(:,:,b,i,1,:) = 1 
            IF( y_pixel(i,1) == sin) bfly(:,:,b,i,2,:) = 1 
         ENDDO
         DO i=1,nbfz 
            IF( z_pixel(i,1) == sin) bflz(:,:,b,i,1,:) = 1 
            IF( z_pixel(i,1) == sin) bflz(:,:,b,i,2,:) = 1 
         ENDDO 
         
         ! 3- solve the boundary-source pixel problem 

          CALL MERGEVOLHCC(nn,npix,nr,nx,ny,nz,
     &                     list_cube_reg,pixel_to_cmpregion,
     &                     aflx(:,:,:,:,oct), Im(1,1,1,1,b,sin,oct))

         DO i=1,nbfx
            sout = x_pixel(i)
            DO bout=1,nb 
               Tm(:,:, bout, sout, b, sin) = 
     &          Tm(:,:, bout, sout, b, sin) + bflx(:,:,bout,i,xout) ! FIXME 
            ENDDO 
         ENDDO
         DO i=1,nbfy
            sout = y_pixel(i)
            DO bout=1,nb 
               Tm(:,:, bout, sout, b, sin) = 
     &          Tm(:,:, bout, sout, b, sin) + bfly(:,:,bout,i,yout) ! FIXME 
            ENDDO 
         ENDDO
         DO i=1,nbfz
            sout = z_pixel(i)
            DO bout=1,nb 
               Tm(:,:, bout, sout, b, sin) = 
     &          Tm(:,:, bout, sout, b, sin) + bflz(:,:,bout,i,zout) ! FIXME 
            ENDDO 
         ENDDO 
      ENDDO 
      ENDDO
         
      END SUBROUTINE HCC_COFHCC_PIXEL


      SUBROUTINE CIRC_CUBE(nr,npix,list_size_reg,list_pix,list_cube_reg)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nr, npix
      INTEGER, INTENT(IN) :: list_size_reg(nr), list_pix(npix,3,nr)

      INTEGER, INTENT(OUT) :: list_cube_reg(2,3,nr)

      INTEGER :: r, maxp


      DO r=1,nr
        maxp = list_size_reg(r)

        list_cube_reg(1,1,r) = MINVAL(list_pix(1:maxp,1,r))
        list_cube_reg(1,2,r) = MINVAL(list_pix(1:maxp,2,r))
        list_cube_reg(1,3,r) = MINVAL(list_pix(1:maxp,3,r))
        list_cube_reg(2,1,r) = MAXVAL(list_pix(1:maxp,1,r))
        list_cube_reg(2,2,r) = MAXVAL(list_pix(1:maxp,2,r))
        list_cube_reg(2,3,r) = MAXVAL(list_pix(1:maxp,3,r))

      END DO

      END SUBROUTINE CIRC_CUBE