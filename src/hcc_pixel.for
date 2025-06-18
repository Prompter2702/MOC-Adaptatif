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
     &                        asrc,
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
      REAL(KIND=8) , DIMENSION(2,nsuo) :: suro_centroid  
      REAL(KIND=8) , DIMENSION(2,nsui) :: suri_centroid 

      REAL(KIND=8) , DIMENSION(3,nr)   :: reg_centroid
      REAL, INTENT(IN)                 :: delt3(3)  ! sides' lengths of the HCC's box 
      
      INTEGER :: pixel_to_cmpregion(npix) ! mapping pixel-to-computational region 
        
      ! outputs 
      REAL(KIND=8),INTENT(INOUT), DIMENSION(ng,ndir,nc,nr,nc,nr)   :: Cm
      REAL(KIND=8),INTENT(INOUT), DIMENSION(ng,ndir,nb,nsuo,nc,nr) :: Em
      REAL(KIND=8),INTENT(INOUT), DIMENSION(ng,ndir,nc,nr,nb,nsui) :: Im
      REAL(KIND=8),INTENT(INOUT), DIMENSION(ng,ndir,nb,nsuo,nb,nsui)::Tm
      
      !! Que des tableaux auxiliaires pour construire les matrices

      REAL, ALLOCATABLE :: asrc(:,:,:,:) ! spatial moments of the angular source 
      REAL, ALLOCATABLE :: aflx(:,:,:,:) ! angular flux along the Z-normal surface 
      REAL, ALLOCATABLE :: bflx(:,:,:,:) ! boundary angular flux along the X-normal surface 
      REAL, ALLOCATABLE :: bfly(:,:,:,:) ! boundary angular flux along the Y-normal surface 
      REAL, ALLOCATABLE :: bflz(:,:,:,:) ! boundary angular flux along the Z-normal surface 
      REAL, ALLOCATABLE :: asrc0(:,:,:)

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
      INTEGER :: normalisation(ndim+1)
      

      REAL  :: asrcreg(ng,ndir,nr,nc)
      
      ! locals 
      INTEGER :: npvol, xin,yin,zin, xout,yout,zout, rin,c,b, rout,
     &           sout, cout, bout, sin, lastadd, nn, i,j,k,p,r

      LOGICAL :: lgki

      REAL(KIND=8), ALLOCATABLE    :: pdslu4(:)
      nn = ng*ndir
      lgki = .FALSE.
      normalisation = (/ 1,nx,ny,nz /)


      ALLOCATE(pdslu4(ndir))
      ALLOCATE(asrc(ng,ndir,npix,nc)) 
      ALLOCATE(aflx(ng,ndir,npix,nc)) 
      ALLOCATE(bflx(ng,ndir,nb,nbfx)) 
      ALLOCATE(bfly(ng,ndir,nb,nbfy)) 
      ALLOCATE(bflz(ng,ndir,nb,nbfz))
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
        asrc0 = 0
        asrc0(:,:,c) = 1.0

        DO p=1,list_size_reg(rin) 
            i = list_pix(p,1,rin)
            j = list_pix(p,2,rin)
            k = list_pix(p,3,rin)
            r = i + (j-1)*nx + (k-1)*nx*ny
            asrc(:,:,r,c) = 1.0/normalisation(c)
        ENDDO

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
     &                   pixel_to_cmpregion,reg_centroid,
     &                   aflx, Cm(1,1,1,1,c,rin))

        CALL MERGEBOUNDHCC(nn,nsuo,nx,ny,nz,
     &                     suro_centroid,   
     &                     nbfx, nbfy, nbfz,
     &                     x_pixel(:,xout), 
     &                     y_pixel(:,yout),
     &                     z_pixel(:,zout),
     &                     bflx(1,1,1,1),
     &                     bfly(1,1,1,1),
     &                     bflz(1,1,1,1),
     &                     Em(1,1,1,1,c,rin))
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

         DO i=1,nbfx 
            IF( x_pixel(i,xin) == sin) bflx(:,:,b,i) = 1/nbfx
         ENDDO
         DO i=1,nbfy 
            IF( y_pixel(i,yin) == sin) bfly(:,:,b,i) = 1/nbfy
         ENDDO
         DO i=1,nbfz 
            IF( z_pixel(i,zin) == sin) bflz(:,:,b,i) = 1/nbfz
         ENDDO 
         
         ! 3- solve the boundary-source pixel problem 
        
         CALL SWEEP_ADA_ONE_OCTANT(nn,ng,ndir,nr,nh,
     &                           nx,ny,nz, delt3,
     &                           xinc,yinc,zinc, oct,
     &                           sigt,mu,eta,ksi,w,pdslu4,
     &                           sgnc,sgni,sgne,sgnt,
     &                           zreg,asrc,aflx,
     &                           bflx, bfly, bflz,tolcor)

          CALL MERGEVOLHCC(nn,npix,nr,nx,ny,nz,
     &                     list_cube_reg,pixel_to_cmpregion,
     &                     aflx, Im(1,1,1,1,b,sin))

          CALL MERGEBOUNDHCC(nn,nsuo,nx,ny,nz,
     &                     suro_centroid,   
     &                     nbfx, nbfy, nbfz,
     &                     x_pixel(:,xout), 
     &                     y_pixel(:,yout),
     &                     z_pixel(:,zout),
     &                     bflx(1,1,1,1),
     &                     bfly(1,1,1,1),
     &                     bflz(1,1,1,1),
     &                     Tm(1,1,1,1,b,sin))
         
      ENDDO 
      ENDDO
         
      END SUBROUTINE HCC_COFHCC_PIXEL



      SUBROUTINE BARY_REG(nr, nx,ny,nz,npix,
     & pixel_to_cmpregion,reg_centroid)

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: nr,npix, nx,ny,nz
        ! mapping pixel-to-computational region
        INTEGER, INTENT(IN) :: pixel_to_cmpregion(npix) 
        REAL(KIND=8), INTENT(INOUT) :: reg_centroid(3,nr) ! centroid of the HCC's regions   

        INTEGER :: r,p,i,j,k

        reg_centroid = 0.0

        DO i=1,nx
        DO j=1,ny
        DO k=1,nz
            p = i + (j-1)*nx + (k-1)*nx*ny
            r= pixel_to_cmpregion(p)
            reg_centroid(1,r) = reg_centroid(1,r) + REAL(i,KIND=8)
            reg_centroid(2,r) = reg_centroid(2,r) + REAL(j,KIND=8)
            reg_centroid(3,r) = reg_centroid(3,r) + REAL(k,KIND=8)
        ENDDO
        ENDDO
        ENDDO

        reg_centroid(1,:)= reg_centroid(1,:)/nx
        reg_centroid(3,:)= reg_centroid(2,:)/ny
        reg_centroid(2,:)= reg_centroid(3,:)/nz

      END SUBROUTINE BARY_REG

      SUBROUTINE BARY_SURF(nsur,nbfx,nbfy,nbfz,
     &                     nx,ny,nz,
     &                     x_pixel,y_pixel,z_pixel,
     &                     xin,yin,zin,nsuo,nsui,
     &                     sur_centroid,suro_centroid,suri_centroid) 


      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nsur,nbfx,nbfy,nbfz
      INTEGER, INTENT(IN) :: nx,ny,nz
      INTEGER, INTENT(IN) :: xin,yin,zin, nsuo,nsui
      INTEGER, DIMENSION(nbfy,2), INTENT(IN) :: y_pixel
      INTEGER, DIMENSION(nbfz,2), INTENT(IN) :: z_pixel
      INTEGER, DIMENSION(nbfy,2), INTENT(IN) :: x_pixel

      REAL(KIND=8), DIMENSION(2,nsur),INTENT(INOUT) :: sur_centroid  
      REAL(KIND=8), DIMENSION(2,nsuo),INTENT(INOUT) :: suro_centroid  
      REAL(KIND=8), DIMENSION(2,nsui),INTENT(INOUT) :: suri_centroid 

      INTEGER :: i,j,k,p

       DO j=1,ny
       DO k=1,nz
        p = j + (k-1)*ny
        suro_centroid(1,x_pixel(p,3-xin)) = 
     &  suro_centroid(1,x_pixel(p,3-xin)) + REAL(i,KIND=8)
     &                      
        suro_centroid(2,x_pixel(p,3-xin)) = 
     &  suro_centroid(2,x_pixel(p,3-xin)) + REAL(j,KIND=8)

        suri_centroid(1,x_pixel(p,xin)) = 
     &  suri_centroid(1,x_pixel(p,xin)) + REAL(i,KIND=8)

        suri_centroid(2,x_pixel(p,xin)) = 
     &  suri_centroid(2,x_pixel(p,xin)) + REAL(j,KIND=8)
       ENDDO
       ENDDO

       DO i=1,nx
       DO k=1,nz
        p = i + (k-1)*ny
        suro_centroid(1,y_pixel(p,3-yin)) = 
     &  suro_centroid(1,y_pixel(p,3-yin)) + REAL(i,KIND=8)
     &                      
        suro_centroid(2,y_pixel(p,3-yin)) = 
     &  suro_centroid(2,y_pixel(p,3-yin)) + REAL(k,KIND=8)

        suri_centroid(1,y_pixel(p,yin)) = 
     &  suri_centroid(1,y_pixel(p,yin)) + REAL(i,KIND=8)

        suri_centroid(2,y_pixel(p,yin)) = 
     &  suri_centroid(2,y_pixel(p,yin)) + REAL(k,KIND=8)
       ENDDO
       ENDDO

       DO i=1,nx
       DO j=1,ny
        p = i + (j-1)*ny
        suro_centroid(1,z_pixel(p,3-zin)) = 
     &  suro_centroid(1,z_pixel(p,3-zin)) + REAL(i,KIND=8)
     &                      
        suro_centroid(2,z_pixel(p,3-zin)) = 
     &  suro_centroid(2,z_pixel(p,3-zin)) + REAL(j,KIND=8)

        suri_centroid(1,z_pixel(p,zin)) = 
     &  suri_centroid(1,z_pixel(p,zin)) + REAL(i,KIND=8)

        suri_centroid(2,z_pixel(p,zin)) = 
     &  suri_centroid(2,z_pixel(p,zin)) + REAL(j,KIND=8)
       ENDDO
       ENDDO

      END SUBROUTINE BARY_SURF 
