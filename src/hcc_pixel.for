       SUBROUTINE HCC_COFHCC_OCT(
    ! inputs HCC main dimensions 
     &                       ng,nreg,nsur,nsui,nsuo,
     &                       nc,nb,nphy,oct,
     &                       x_pixel,y_pixel,z_pixel,
     &                       sur_centroid,
     &                       reg_centroid,
     &                       pixel_to_cmpregion, iosurf_map, 
     &                       list_size_reg,cmpreg_to_pixel,size_sum_reg,
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
     &                       asrc,
       ! input: pixel-to-medium array 
     &                       zreg, 
       ! output: angular flux moments 
     &                       aflx,
       ! Delta on each direction of the region
     &                       delt3,
      ! max inner iterations and tolerance for inner iterations
     &                       tolcor,
      ! outputs  
     &                       Cm,Im,Em,Tm,
     &                       pdslu4,
     &                       asrc0,
     &                       aflx0, aflx1,
     &                       asrcm0, asrcm1,
     &                       finc1,fout1,
     &                       finc0,fout0,
     &                       ccof,icof,
     &                       ecof,tcof,
     &                       ccof8,icof8,
     &                       ecof8,tcof8,
     &                       aflxmean,
     &                       xstlvl1,
     &                       errcor,errmul)

      USE FLGINC
      USE SRCCORR
      USE SWEEP8ONE
      
      IMPLICIT NONE

      INTEGER, PARAMETER :: n8 = 8, ns4 = 4, ndim = 3
      
      ! HCC data inputs 
      INTEGER :: ng      ! n. of groups 
      INTEGER :: nreg      ! n. of regions in HCC
      INTEGER :: nsur    ! n. of surfaces of the HCC
      INTEGER :: nsui    ! n. of incoming surfaces of the HCC
      INTEGER :: nsuo    ! n. of outgoing surfaces of the HCC
      INTEGER :: nc      ! n. of spatial moments per HCC's volume
      INTEGER :: nb      ! n. of spatial moments per HCC's surface 
      INTEGER :: nphy    ! n. of physical HCCs sharing the same HCC geometry

      ! data for pixels in the HCC (local data)
      
      INTEGER :: npix    ! n. of pixels -> nreg 
      INTEGER :: nbd     ! nb * (ndim = 3) 
      INTEGER :: nmat    ! n. of media 
      INTEGER :: nbfx    ! n. of pixel surfaces in X direction
      INTEGER :: nbfy    ! n. of pixel surfaces in Y direction
      INTEGER :: nbfz    ! n. of pixel surfaces in Z direction
      INTEGER :: nx      ! n. of pixels in X direction
      INTEGER :: ny      ! n. of pixels in Y direction
      INTEGER :: nz      ! n. of pixels in Z direction

      INTEGER :: list_size_reg(nreg)  ! list of nb of pixel in each region
      INTEGER :: list_cube_reg(2,3,nreg) ! list of coordinates of the circonsrit 
      ! cube of each region
      INTEGER, INTENT(IN) :: x_pixel(nbfx,2) ! list of pixel surfaces in X direction
      INTEGER, INTENT(IN) :: y_pixel(nbfy,2) ! list of pixel surfaces in Y direction
      INTEGER, INTENT(IN) :: z_pixel(nbfz,2) ! list of pixel surfaces in Z direction

      INTEGER :: list_oct_surf(8,nsur) ! list of surfaces per incomming octant

      INTEGER, INTENT(IN) :: oct      ! octant index
      REAL    , DIMENSION(ng,nmat)    :: sigt   ! total XS per HCC's region
      REAL    , DIMENSION(ng,nphy,nreg) :: sigs   ! total XS per HCC's region
      
      ! incoming/outgoing surface index per HCC's surface (for the octant)
      ! First the incoming surfaces, then the outgoing surfaces
      INTEGER, INTENT(IN), DIMENSION(nsur,2) :: iosurf_map
       
      REAL(KIND=8) , DIMENSION(2,nsur) :: sur_centroid  
      REAL(KIND=8) , DIMENSION(2,nsur) :: suro_centroid  
      REAL(KIND=8) , DIMENSION(2,nsur) :: suri_centroid  

      REAL(KIND=8) , DIMENSION(3,nreg)   :: reg_centroid
      REAL, INTENT(IN)                 :: delt3(3)  ! sides' lengths of the HCC's box 
      
      INTEGER :: pixel_to_cmpregion(npix) ! mapping pixel-to-computational region 
        
      ! outputs 
      REAL(KIND=8),INTENT(INOUT),
     &                 DIMENSION(ng,ndir,nc,nreg,nc,nreg)  :: Cm
      REAL(KIND=8),INTENT(INOUT),
     &                 DIMENSION(ng,ndir,nb,nsuo,nc,nreg)  :: Em
      REAL(KIND=8),INTENT(INOUT),
     &                 DIMENSION(ng,ndir,nc,nreg,nb,nsui)  :: Im
      REAL(KIND=8),INTENT(INOUT),
     &                 DIMENSION(ng,ndir,nb,nsuo,nb,nsui)  ::Tm
      
      !! Que des tableaux auxiliaires pour construire les matrices


      REAL(KIND=8), INTENT(INOUT) :: pdslu4(ndir)
      REAL, INTENT(INOUT) :: asrc(ng,ndir,npix,nc)
      REAL, INTENT(INOUT) :: aflx(ng,ndir,npix,nc)
      REAL, INTENT(INOUT) :: bflx(ng,ndir,nb,nbfx)
      REAL, INTENT(INOUT) :: bfly(ng,ndir,nb,nbfy)
      REAL, INTENT(INOUT) :: bflz(ng,ndir,nb,nbfz)
      REAL, INTENT(INOUT) :: asrc0(ng,ndir,nc)


      REAL, INTENT(INOUT) :: aflx0(ng*ndir,nc), aflx1(ng*ndir,nc,n8)
      REAL, INTENT(INOUT) :: asrcm0(ng,ndir,nc), asrcm1(ng,ndir,nc,n8)
      REAL, INTENT(INOUT) :: finc1(ng*ndir,nb,3,ns4)
      REAL, INTENT(INOUT) :: fout1(ng*ndir,nb,3,ns4)
      REAL, INTENT(INOUT) :: finc0(ng*ndir,nb,3), fout0(ng*ndir,nb,3)
      REAL, INTENT(INOUT) :: ccof(ng*ndir,nc,nc),icof(ng*ndir,nc,nbd)
      REAL, INTENT(INOUT) :: ecof(ng*ndir,nbd,nc),tcof(ng*ndir,nbd,nbd)
      REAL, INTENT(INOUT) :: ccof8(ng*ndir,nc,nc,n8),
     &                       icof8(ng*ndir,nc,nbd,n8),
     &                       ecof8(ng*ndir,nbd,nc,n8),
     &                       tcof8(ng*ndir,nbd,nbd,n8)
      REAL, INTENT(INOUT) :: aflxmean(ng*ndir,npix)
      REAL, INTENT(INOUT) :: xstlvl1(ng,n8)
      REAL(KIND=8), INTENT(INOUT) :: errcor(ng, ndir),errmul(ng, ndir)
      INTEGER, INTENT(IN) :: cmpreg_to_pixel(npix),size_sum_reg(nreg)


      INTEGER, INTENT(IN) :: zreg(npix) ! pixel-to-medium assignment      
      
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
      REAL, INTENT(IN)         :: sphr(nhrm,nd)
      
      !  sign matrices 
      INTEGER,INTENT(IN) :: sgnc(nc,nc),sgni(nc,nbd)
      INTEGER,INTENT(IN) :: sgne(nbd,nc),sgnt(nbd,nbd)
      REAL, INTENT(IN) :: tolcor
      INTEGER :: normalisation(ndim+1)
      

      REAL  :: asrcreg(ng,ndir,nreg,nc), bflxsur(ng,ndir,nb)
      INTEGER :: nb_cell
      LOGICAL :: ok,oksrc
      REAL :: errtot,errbnd(3)
      ! locals 
      INTEGER :: npvol, xin,yin,zin, xout,yout,zout, rin,c,b, rout,
     &           sout, cout, bout, sin, lastadd, nn, p, cnt,g, i,j,k,r
      LOGICAL :: lgki


      nn = ng*ndir
      lgki = .FALSE.
      normalisation = (/ 1,nx,ny,nz /)

      nb_cell = 0
      ok =.FALSE.
      oksrc = .FALSE.
      errtot = 0.0
      errbnd = 0.0

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      xin = xinc(oct)
      yin = yinc(oct)
      zin = zinc(oct)
      xout = 3-xin
      yout = 3-yin 
      zout = 3-zin 

      Cm = 0.0
      Em = 0.0
      Im = 0.0
      Tm = 0.0

      ! n. of volume-source problems 
      
      DO rin=1,nreg
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

        CALL SPLIT_SRC_REG(ng*ndir,npix,nx,ny,nz,list_size_reg(rin),rin,
     &                     asrc0,asrc,cmpreg_to_pixel,size_sum_reg(rin),
     &                     reg_centroid(:,rin))

    ! 3- solve the pixel volume source problem 
    ! Construire une fonction plus adaptée !!
        ok =.FALSE.
        oksrc = .FALSE.
        errtot = 0.0
        errbnd = 0.0

        CALL SRCHOMO0(ng*ndir,npix,nx,ny,
     &                 1,nx,1,ny,1,nz,
     &                 zreg,asrc,asrcm0)
        CALL ONE_COEF3D(ng*ndir,ndir,ng,mu,eta,ksi,
     &                delt3,sigt(:,zreg(1)),
     &                sgnc,sgni,sgne,sgnt,
     &                ccof,icof,ecof,tcof)
        CALL MERGEBOUND0(ng*ndir,nb,
     &                   1,nx,1,ny,1,nz,
     &                   nx,ny,nz,
     &                   bflx, bfly,bflz,finc0)

        CALL SWEEP_ONEREGION(ng*ndir,2,asrcm0, finc0, aflx0, fout0,
     &                          ccof,icof,ecof,tcof)

       CALL RECADA_ONE_OCTANT(ng*ndir,ng,ndir,npix,nh,
     &                        nx,ny,nz,
     &                        xinc,yinc,zinc, oct,
     &                        1,nx,1,ny,1,nz,0,
     &                        sigt, xstlvl1,
     &                        mu,eta,ksi,w,pdslu4,
     &                        sgnc,sgni,sgne,sgnt,
     &                        zreg,asrc,aflx,
     &                        bflx, bfly, bflz,
     &                        aflx0, aflx1,
     &                        asrcm0, asrcm1,finc0,finc1,
     &                        fout0,fout1,ccof,icof,ecof,tcof,
     &                        ccof8,icof8,ecof8,tcof8,tolcor,
     &                        errcor,errmul,ok,delt3,i,j,k,r,
     &                        nb_cell,g,cnt,oksrc,errtot,aflxmean,
     &                        errbnd)

        ! 4- fill-in the C and E matrix column
        CALL MERGEVOLHCC(ng,ndir,npix,nreg,nx,ny,nz,
     &                   pixel_to_cmpregion,reg_centroid,
     &                   aflx, Cm(1,1,1,1,c,rin) )
        
        CALL MERGEBOUNDHCC(nn,nsui,nsuo,nx,ny,nz,
     &                     sur_centroid,iosurf_map,
     &                     nbfx, nbfy, nbfz,
     &                     x_pixel(:,xout), 
     &                     y_pixel(:,yout),
     &                     z_pixel(:,zout),
     &                     bflx, bfly,bflz,
     &                     Em(1,1,1,1,c,rin))
        END DO


      ENDDO 
          
      ! n. of surface-source problems 

      DO sin=1,nsur
      IF (iosurf_map(sin,2) == 2) CYCLE ! skip if no incoming surface
      
      DO b=1,nb 

         ! 1- contruct boundary-source for the couple of indexes (s,b) 
         ! s = HCC surface, b = spatial surface component 
         ! 2- projection: project the boundary HCC source onto
         ! the pixel boundary support  
        asrc = 0.0
        bflx = 0.0
        bfly = 0.0
        bflz = 0.0
        asrc0 = 0.0
        bflxsur = 0.0
        bflxsur(:,:,b) = 1.0

        CALL SPLIT_BOUND_SUR(nn,nx,ny,nz,sin,
     &                       x_pixel(:,xin),
     &                       y_pixel(:,yin),
     &                       z_pixel(:,zin),bflxsur,
     &                       bflx,bfly,bflz, sur_centroid(:,sin))
         
        ! 3- solve the boundary-source pixel problem 
        ok =.FALSE.
        oksrc = .FALSE.
        errtot = 0.0
        errbnd = 0.0

        CALL SRCHOMO0(nn,npix,nx,ny,
     &                 1,nx,1,ny,1,nz,
     &                 zreg,asrc,asrcm0)
        CALL ONE_COEF3D(nn,ndir,ng,mu,eta,ksi,
     &                delt3,sigt(:,zreg(1)),
     &                sgnc,sgni,sgne,sgnt,
     &                ccof,icof,ecof,tcof)
        CALL MERGEBOUND0(nn,nb,
     &                   1,nx,1,ny,1,nz,
     &                   nx,ny,nz,
     &                   bflx, bfly,bflz,finc0)
        CALL SWEEP_ONEREGION(nn,2,asrcm0, 
     &                         finc0, aflx0, fout0,
     &                         ccof,icof,ecof,tcof)

        CALL RECADA_ONE_OCTANT(nn,ng,ndir,npix,nh,
     &                       nx,ny,nz,
     &                       xinc,yinc,zinc, oct,
     &                       1,nx,1,ny,1,nz,0,
     &                       sigt, xstlvl1,
     &                       mu,eta,ksi,w,pdslu4,
     &                       sgnc,sgni,sgne,sgnt,
     &                       zreg,asrc,aflx,
     &                       bflx, bfly, bflz,
     &                       aflx0, aflx1,
     &                       asrcm0, asrcm1,finc0,finc1,
     &                       fout0,fout1,ccof,icof,ecof,tcof,
     &                       ccof8,icof8,ecof8,tcof8,tolcor,
     &                       errcor,errmul,ok,delt3,i,j,k,r,
     &                       nb_cell,g,cnt,oksrc,errtot,aflxmean,
     &                       errbnd)

          CALL MERGEVOLHCC(ng,ndir,npix,nreg,nx,ny,nz,
     &                     pixel_to_cmpregion,reg_centroid,
     &                     aflx, Im(1,1,1,1,b,iosurf_map(sin,1)))


          CALL MERGEBOUNDHCC(nn,nsuo,nsui,nx,ny,nz,
     &                     sur_centroid, iosurf_map,
     &                     nbfx, nbfy, nbfz,
     &                     x_pixel(:,xout), 
     &                     y_pixel(:,yout),
     &                     z_pixel(:,zout),
     &                     bflx(1,1,1,1),
     &                     bfly(1,1,1,1),
     &                     bflz(1,1,1,1),
     &                     Tm(1,1,1,1,b,iosurf_map(sin,1)))

         
      ENDDO 
      ENDDO
         
      END SUBROUTINE HCC_COFHCC_OCT


      !-----------------------------------------------------------------

      
      SUBROUTINE HCC_COEFF(ng,nn,nreg,nsur,nsui,nsuo,
     &                     nc,nb,nphy, x_pixel,y_pixel,z_pixel,
     &                     pixel_to_cmpregion,iosurf_map,
     &                     list_size_reg,
     &                     npix,nh,nmat,
     &                     nbd,nbfx,nbfy,nbfz,
     &                     nx,ny,nz,
     &                     nani,nhrm,nd,ndir,noct,
     &                     sigt,sigs,
     &                     asrc,aflx,bflx, bfly, bflz, 
     &                     mu,eta,ksi,w,pisn,sphr,
     &                     zreg, sgnc,sgni, sgne,sgnt,
     &                     delt3, pdslu4,aflxmean,
     &                     tolcor,asrc0,
     &                     aflx0, aflx1,
     &                     asrcm0, asrcm1,
     &                     finc1, fout1,
     &                     finc0, fout0,
     &                     ccof,icof,
     &                     ecof,tcof,
     &                     ccof8,icof8,
     &                     ecof8,tcof8,Cm,Em,Im,Tm,errcor,errmul)

      IMPLICIT NONE

      INTEGER, PARAMETER :: n8=8, ns4=4

      INTEGER, INTENT(IN) :: ng,nreg,nsur,nsui,nsuo,nn
      INTEGER, INTENT(IN) :: nc,nb,nphy
      INTEGER, INTENT(IN) :: pixel_to_cmpregion(npix) ! mapping pixel-to-computational region
      INTEGER, INTENT(IN) :: npix,nh,nmat
      INTEGER, INTENT(IN) :: nbd,nbfx,nbfy,nbfz
      INTEGER, INTENT(IN) :: nx,ny,nz
      INTEGER, INTENT(IN) :: nani,nhrm,nd,ndir,noct
      
      REAL, INTENT(IN) :: sigt(ng, nmat), sigs(ng,0:nani,nmat)
      REAL(KIND=8), INTENT(IN) :: mu(ndir),eta(ndir),ksi(ndir),
     &                            w(ndir),pisn
      REAL, INTENT(IN)  :: sphr(nhrm,nd)
      INTEGER, INTENT(IN) ::  zreg(npix) 
      REAL, INTENT(IN) :: delt3(3),tolcor
      INTEGER, INTENT(IN) :: iosurf_map(nsur,2,noct)
      INTEGER, INTENT(IN) :: sgnc(nc,nc,noct),sgni(nc,nbd,noct)
      INTEGER, INTENT(IN) :: sgne(nc,nbd,noct),sgnt(nbd,nbd,noct)
      
      
      INTEGER, INTENT(INOUT) :: list_size_reg(nreg)
      REAL, INTENT(INOUT) :: asrc(ng,ndir,npix,nc) 
      REAL, INTENT(INOUT) :: aflx(ng,ndir,npix,nc) 
      REAL, INTENT(INOUT) :: bflx(ng,ndir,nb,nbfx) 
      REAL, INTENT(INOUT) :: bfly(ng,ndir,nb,nbfy) 
      REAL, INTENT(INOUT) :: bflz(ng,ndir,nb,nbfz)

      
      REAL, INTENT(INOUT)   :: asrc0(ng,ndir,nc)
      REAL, INTENT(INOUT) :: aflx0(nn,nc), aflx1(nn,nc,n8)
      REAL, INTENT(INOUT) :: asrcm0(ng,ndir,nc), asrcm1(ng,ndir,nc,n8)
      REAL, INTENT(INOUT) :: finc1(nn,nb,3,ns4), fout1(nn,nb,3,ns4)
      REAL, INTENT(INOUT) :: finc0(nn,nb,3),    fout0(nn,nb,3)
      REAL, INTENT(INOUT) :: ccof(nn,nc,nc),icof(nn,nc,nbd)
      REAL, INTENT(INOUT) :: ecof(nn,nbd,nc),tcof(nn,nbd,nbd)
      REAL, INTENT(INOUT) :: ccof8(nn,nc,nc,n8),icof8(nn,nc,nbd,n8)
      REAL, INTENT(INOUT) :: ecof8(nn,nbd,nc,n8),tcof8(nn,nbd,nbd,n8)


      REAL, INTENT(INOUT) :: aflxmean(nn,npix,noct)
      REAL(KIND=8), INTENT(INOUT)  :: errcor(ng, ndir),errmul(ng, ndir)
      

      REAL(KIND=8), INTENT(INOUT)  :: Cm(ng,ndir,nc, nreg,nc,nreg,noct) 
      REAL(KIND=8), INTENT(INOUT)  :: Em(ng,ndir,nb, nsuo,nc,nreg,noct) 
      REAL(KIND=8), INTENT(INOUT)  :: Im(ng,ndir,nc, nreg,nb,nsui,noct) 
      REAL(KIND=8), INTENT(INOUT)  :: Tm(ng,ndir,nb, nsuo,nb,nsui,noct)


      REAL(KIND=8), DIMENSION(3,nreg) :: reg_centroid 
      REAL(KIND=8), DIMENSION(2,nsur) :: sur_centroid  
      REAL, DIMENSION(ng,n8)  :: xstlvl1


      INTEGER,INTENT(IN) :: x_pixel(nbfx,2) ! list of pixel surfaces in X direction
      INTEGER,INTENT(IN) :: y_pixel(nbfy,2) ! list of pixel surfaces in Y direction
      INTEGER,INTENT(IN) :: z_pixel(nbfz,2) ! list of pixel surfaces in Z direction
      REAL(KIND=8), INTENT(INOUT) :: pdslu4(ndir)

      INTEGER :: oct
      INTEGER, ALLOCATABLE :: cmpreg_to_pixel(:), size_sum_reg(:) ! mapping computational region to pixel

      
      ALLOCATE(cmpreg_to_pixel(npix))
      ALLOCATE(size_sum_reg(nreg))


      CALL BARY_REG(nreg, nx,ny,nz,npix,pixel_to_cmpregion,
     &              reg_centroid,list_size_reg)

      CALL CMPREG_TO_PIX(nreg,npix, list_size_reg, pixel_to_cmpregion, 
     &                        size_sum_reg,cmpreg_to_pixel )

      CALL BARY_SURF(nsur,nbfx,nbfy,nbfz,
     &                     nx,ny,nz,
     &                     x_pixel,y_pixel,z_pixel,
     &                     sur_centroid) 

      DO oct=1,noct

         CALL HCC_COFHCC_OCT(ng,nreg,nsur,nsui,nsuo,
     &                       nc,nb,nphy,oct,
     &                       x_pixel,y_pixel,z_pixel,
     &                       sur_centroid,
     &                       reg_centroid,
     &                       pixel_to_cmpregion, iosurf_map(1,1,oct),
     &                       list_size_reg,cmpreg_to_pixel,size_sum_reg,
     &                       npix,nh,nmat,
     &                       nbd,nbfx,nbfy,nbfz,
     &                       nx,ny,nz,
     &                       nani,nhrm,nd,ndir,noct,
     &                       sigt,sigs,
     &                       bflx,bfly,bflz,
     &                       mu,eta,ksi,w,pisn,sphr,
     &                       sgnc(:,:,oct),sgni(:,:,oct),
     &                       sgne(:,:,oct),sgnt(:,:,oct),
     &                       asrc,zreg, aflx,delt3,tolcor,
     &                       Cm(:,:,:, :,:,:,oct),
     &                       Im(:,:,:, :,:,:,oct),
     &                       Em(:,:,:, :,:,:,oct),
     &                       Tm(:,:,:, :,:,:,oct),
     &                       pdslu4,
     &                       asrc0,aflx0, aflx1,
     &                       asrcm0, asrcm1,
     &                       finc1,fout1,finc0,fout0,
     &                       ccof,icof,ecof,tcof,
     &                       ccof8,icof8,ecof8,tcof8,
     &                       aflxmean(1,1,oct),
     &                       xstlvl1,
     &                       errcor,errmul)

      ENDDO


      END SUBROUTINE HCC_COEFF

      !-----------------------------------------------------------------


      SUBROUTINE BARY_REG(nreg, nx,ny,nz,npix,pixel_to_cmpregion,
     &                    reg_centroid,list_size_reg)

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: nreg,npix, nx,ny,nz
        ! mapping pixel-to-computational region
        INTEGER, INTENT(IN) :: pixel_to_cmpregion(npix) 
        REAL(KIND=8), INTENT(INOUT) :: reg_centroid(3,nreg) ! centroid of the HCC's regions   
        INTEGER, INTENT(INOUT) :: list_size_reg(nreg) ! list of nb of pixel in each region

        INTEGER :: r,p,i,j,k

        reg_centroid = 0.0
        list_size_reg = 0

        DO i=1,nx
        DO j=1,ny
        DO k=1,nz
            p = i + (j-1)*nx + (k-1)*nx*ny
            r= pixel_to_cmpregion(p)
            list_size_reg(r) = list_size_reg(r) + 1
            reg_centroid(1,r) = reg_centroid(1,r) + REAL(i,KIND=8)
            reg_centroid(2,r) = reg_centroid(2,r) + REAL(j,KIND=8)
            reg_centroid(3,r) = reg_centroid(3,r) + REAL(k,KIND=8)
        ENDDO
        ENDDO
        ENDDO

        reg_centroid(1,:)= reg_centroid(1,:)/list_size_reg(:)
        reg_centroid(2,:)= reg_centroid(2,:)/list_size_reg(:)
        reg_centroid(3,:)= reg_centroid(3,:)/list_size_reg(:)

      END SUBROUTINE BARY_REG


      SUBROUTINE CMPREG_TO_PIX(nreg,npix, list_size_reg, 
     &                         pixel_to_cmpregion,
     &                         size_sum, cmpreg_to_pixel )

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nreg,npix, list_size_reg(nreg),
     &                       pixel_to_cmpregion(npix)
      INTEGER, INTENT(OUT) :: cmpreg_to_pixel(npix)

      INTEGER :: r,p, cnt(nreg), size_sum(nreg)

      size_sum = 0
      cnt = 0

      DO r=2,nreg
         size_sum(r) = size_sum(r-1) + list_size_reg(r-1)
      ENDDO

      DO p=1,npix
         r = pixel_to_cmpregion(p)
         cnt(r) = cnt(r) + 1
         cmpreg_to_pixel(cnt(r)+size_sum(r)) = p
      END DO
      

      END SUBROUTINE

      SUBROUTINE BARY_SURF(nsur,nbfx,nbfy,nbfz,
     &                     nx,ny,nz,
     &                     x_pixel,y_pixel,z_pixel,
     &                     sur_centroid) 


      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nsur,nbfx,nbfy,nbfz
      INTEGER, INTENT(IN) :: nx,ny,nz
      INTEGER, DIMENSION(nbfx,2), INTENT(IN) :: x_pixel
      INTEGER, DIMENSION(nbfy,2), INTENT(IN) :: y_pixel
      INTEGER, DIMENSION(nbfz,2), INTENT(IN) :: z_pixel
      
      REAL(KIND=8), DIMENSION(2,nsur),INTENT(INOUT) :: sur_centroid  

      INTEGER :: i,j,k,p, cnt(nsur)
      
      cnt = 0
      sur_centroid = 0.0

      DO j=1,ny
      DO k=1,nz
        p = j + (k-1)*ny
        cnt(x_pixel(p,1)) = cnt(x_pixel(p,1)) + 1
        cnt(x_pixel(p,2)) = cnt(x_pixel(p,2)) + 1

        sur_centroid(1,x_pixel(p,1)) = 
     &  sur_centroid(1,x_pixel(p,1)) + REAL(j,KIND=8)
        sur_centroid(2,x_pixel(p,1)) = 
     &  sur_centroid(2,x_pixel(p,1)) + REAL(k,KIND=8)
        
        sur_centroid(1,x_pixel(p,2)) = 
     &  sur_centroid(1,x_pixel(p,2)) + REAL(j,KIND=8)
        sur_centroid(2,x_pixel(p,2)) = 
     &  sur_centroid(2,x_pixel(p,2)) + REAL(k,KIND=8)
      ENDDO
      ENDDO

      DO i=1,nx
      DO k=1,nz
        p = i + (k-1)*ny
        cnt(y_pixel(p,1)) = cnt(y_pixel(p,1)) + 1
        cnt(y_pixel(p,2)) = cnt(y_pixel(p,2)) + 1
        
        sur_centroid(1,y_pixel(p,1)) = 
     &  sur_centroid(1,y_pixel(p,1)) + REAL(i,KIND=8)
        sur_centroid(2,y_pixel(p,1)) = 
     &  sur_centroid(2,y_pixel(p,1)) + REAL(k,KIND=8)

        sur_centroid(1,y_pixel(p,2)) = 
     &  sur_centroid(1,y_pixel(p,2)) + REAL(i,KIND=8)
        sur_centroid(2,y_pixel(p,2)) = 
     &  sur_centroid(2,y_pixel(p,2)) + REAL(k,KIND=8)
      ENDDO
      ENDDO

      DO i=1,nx
      DO j=1,ny
        p = i + (j-1)*ny
        cnt(z_pixel(p,1)) = cnt(z_pixel(p,1)) + 1
        cnt(z_pixel(p,2)) = cnt(z_pixel(p,2)) + 1

        sur_centroid(1,z_pixel(p,1)) = 
     &  sur_centroid(1,z_pixel(p,1)) + REAL(i,KIND=8)
        sur_centroid(2,z_pixel(p,1)) = 
     &  sur_centroid(2,z_pixel(p,1)) + REAL(j,KIND=8)

        sur_centroid(1,z_pixel(p,2)) = 
     &  sur_centroid(1,z_pixel(p,2)) + REAL(i,KIND=8)
        sur_centroid(2,z_pixel(p,2)) = 
     &  sur_centroid(2,z_pixel(p,2)) + REAL(j,KIND=8)
      ENDDO
      ENDDO

      sur_centroid(1,:) = sur_centroid(1,:) / cnt(:)
      sur_centroid(2,:) = sur_centroid(2,:) / cnt(:)

      END SUBROUTINE BARY_SURF 
