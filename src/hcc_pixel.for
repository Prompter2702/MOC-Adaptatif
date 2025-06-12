       SUBROUTINE HCC_COFHCC_PIXEL(
     &                       oct,
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
       ! source moments and self-scattering xs 
     &                       sigg, !tmom,
       ! input: pixel-to-medium array 
     &                       zreg, 
       ! output: angular flux moments 
     &                       aflx,
       ! Delta on each direction of the region
     &                       delt3,
      ! max inner iterations and tolerance for inner iterations
     &                       maxinner, tolinner,
      ! outputs  
     &                       Cm,Im,Em,Tm )

      USE FLGINC
      
      IMPLICIT NONE
      
      ! HCC data inputs 
      INTEGER :: ng      ! n. of groups 
      INTEGER :: nr      ! n. of regions in the HCC 
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
      INTEGER :: list_pix(npix,nr)  ! list of pixels in each region 
      INTEGER :: x_pixel(nbfx) ! list of pixel surfaces in X direction
      INTEGER :: y_pixel(nbfy) ! list of pixel surfaces in Y direction
      INTEGER :: z_pixel(nbfz) ! list of pixel surfaces in Z direction 

      INTEGER :: oct      ! octant index
      REAL    , DIMENSION(ng,nmat)    :: sigt   ! total XS per HCC's region
      REAL    , DIMENSION(ng,nphy,nr) :: sigs   ! total XS per HCC's region
     
      INTEGER , DIMENSION(nsur) :: iosurf_map  ! incoming/outgoing surface index per HCC's surface (per octant )
       
      REAL(KIND=8) , DIMENSION(2,nsur) :: sur_centroid  
      REAL(KIND=8) , DIMENSION(3,nr)   :: reg_centroid
      REAL, INTENT(IN)                 :: delt3(3)  ! sides' lengths of the HCC's box 
      
      INTEGER :: pixel_to_cmpregion(npix) ! mapping pixel-to-computational region 
        
      ! outputs 
      REAL(KIND=8), DIMENSION(ng*nphy,ndir, nc, nr  ,nc,nr) :: Cm
      REAL(KIND=8), DIMENSION(ng*nphy,ndir, nb, nsuo,nc,nr) :: Em
      
      REAL(KIND=8), DIMENSION(ng*nphy,ndir, nc, nr  , nb, nsui) :: Im
      REAL(KIND=8), DIMENSION(ng*nphy,ndir, nb, nsuo, nb, nsui) :: Tm
      
      
      REAL :: flxm(ng*nphy,npix,nh,nc) ! angular moments of the flux 
      REAL :: asrc(ng*nphy,ndir,npix,nc) ! spatial moments of the angular source 
      REAL :: aflx(ng*nphy,ndir,npix,nc) ! angular flux along the Z-normal surface 
      
      REAL :: bflx(ng*nphy,ndir,nb,nbfx,2)  ! boundary angular flux along the X-normal surface 
      REAL :: bfly(ng*nphy,ndir,nb,nbfy,2)  ! boundary angular flux along the Y-normal surface 
      REAL :: bflz(ng*nphy,ndir,nb,nbfz,2)  ! boundary angular flux along the Z-normal surface 
      INTEGER, INTENT(IN) :: zreg(nphy,npix)     ! pixel-to-medium assignment      
      
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
      INTEGER , INTENT(IN) :: maxinner
      REAL, INTENT(IN) :: tolinner


      REAL  :: asrcreg(ng*nphy,ndir, nc, nr)
      
      ! locals 
      INTEGER :: npvol, xin,yin,zin, xout,yout,zout, rin,c,b, rout,i,
     &           sout, cout, bout, sin
      
      !  ngp = nphy * ng 
      !  ngd = ndir * ng    
      !  ngdp = nphy * ngd
      
      xin = xinc(oct)
      xout = 3-xin
      yin = yinc(oct)
      yout = 3-yin 
      zin = zinc(oct)
      zout = 3-zin 

      ! n. of volume-source problems 
      DO rin=1,nr
      DO c=1,nc 
         
        ! 1- contruct the volume-source for the couple of indexes (r,c) r = HCC region, c = spatial volume component 
        ! 2- projection: project the volume HCC source onto the pixel support  
        asrc = 0
        bflx = 0
        bfly = 0
        bflz = 0


        CALL SPLITVOLHCC(ng*nphy*ndir,npix, list_size_reg(rin),
     &                  nx,ny,nz,
     &                  1,nx,1,ny,1,nz,
     &                  asrcreg(:,:,:,rin) ,asrc, reg_centroid(:,rin),
     &                  list_pix(:,rin))

        ! 3- solve the pixel volume source problem 
    ! Construire une fonction plus adaptée !!
    !     CALL SWEE3D_ADAPTIVE(ng*ndir,ng,nr,nh,nc,nmat,
    !  &                       nb,nbd,nbfx,nbfy,nbfz,
    !  &                       nx,ny,nz,nani,nhrm,nd,ndir,
    !  &                       sigt,sigs,srcm,
    !  &                       bflx,bfly,bflz,
    !  &                       mu,eta,ksi,w,pisn,sphr,
    !  &                       sgnc,sgni,sgne,sgnt,sigg,
    !  &                       dsrc,rdir,zreg, 
    !  &                       dira,dirf,lgki,
    !  &                       flxm,delt3,
    !  &                       maxinner, tolinner)
        
        ! 4- fill-in the C and E matrix column
        DO i=1,npix
            rout = pixel_to_cmpregion(i)
            DO cout=1,nc 
                Cm(:, :, cout, rout,c,rin) = 
     &          Cm(:, :, cout, rout,c,rin) + aflx(:,:,i,c) ! FIXME 
            ENDDO  
        ENDDO 

        DO i=1,nbfx
           sout = x_pixel(i)
           DO bout=1,nb 
              Em(:, :, bout, sout,c,rin) = 
     &        Em(:, :, bout, sout,c,rin) + bflx(:,:,bout,i,xout) ! FIXME 
           ENDDO 
        ENDDO 

        DO i=1,nbfy
           sout = y_pixel(i)
           DO bout=1,nb 
              Em(:, :, bout, sout,c,rin) = 
     &        Em(:, :, bout, sout,c,rin) + bfly(:,:,bout,i,yout) ! FIXME 
           ENDDO 
        ENDDO 

        DO i=1,nbfz
           sout = z_pixel(i)
           DO bout=1,nb 
              Em(:, :, bout, sout,c,rin) = 
     &        Em(:, :, bout, sout,c,rin) + bflz(:,:,bout,i,zout) ! FIXME 
           ENDDO 
        ENDDO 
        
      ENDDO 
      ENDDO 
          
      ! n. of surface-source problems 
      DO sin=1,nsui
      DO b=1,nb 
         ! 1- contruct boundary-source for the couple of indexes (s,b) 
         ! s = HCC surface, b = spatial surface component 
      
         ! 2- projection: project the boundary HCC source onto
         ! the pixel boundary support  
         asrc = 0
         bflx = 0
         bfly = 0
         bflz = 0    
        
         DO i=1,npix 
            IF( x_pixel(i) == sin) bflx(:,:,b,i,xinc) = 1 ! FIXME 
         ENDDO
         DO i=1,npix 
            IF( y_pixel(i) == sin) bfly(:,:,b,i,yinc) = 1 ! FIXME 
         ENDDO
         DO i=1,npix 
            IF( z_pixel(i) == sin) bflz(:,:,b,i,zinc) = 1 ! FIXME 
         ENDDO 
         
         ! 3- solve the boundary-source pixel problem 
        !  CALL SWEE3D_ADAPTIVE
         
         ! 4- fill_in the I and T matrix colums
         DO i=1,npix
             rout = pixel_to_cmpregion(i)
             DO cout=1,nc 
                Im(:,:, cout, rout, b, sin) = 
     &            Im(:,:, cout, rout, b, sin) + aflx(:,:,i,c) ! FIXME 
             ENDDO  
         ENDDO 
         DO i=1,nbfx
            sout = x_pixel(i)
            DO bout=1,nb 
               Tm(:,:, bout, sout, b, sin) = 
     &         Tm(:,:, bout, sout, b, sin) + bflx(:,:,bout,i,xout) ! FIXME 
            ENDDO 
         ENDDO
         DO i=1,nbfy
            sout = y_pixel(i)
            DO bout=1,nb 
               Tm(:,:, bout, sout, b, sin) = 
     &         Tm(:,:, bout, sout, b, sin) + bfly(:,:,bout,i,yout) ! FIXME 
            ENDDO 
         ENDDO
         DO i=1,nbfz
            sout = z_pixel(i)
            DO bout=1,nb 
               Tm(:,:, bout, sout, b, sin) = 
     &         Tm(:,:, bout, sout, b, sin) + bflz(:,:,bout,i,zout) ! FIXME 
            ENDDO 
         ENDDO 
      ENDDO 
      ENDDO
         
      END SUBROUTINE HCC_COFHCC_PIXEL