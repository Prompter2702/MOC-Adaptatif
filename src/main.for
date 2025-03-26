!****************************************************************************
!     "nb"    is nubmer of spatial moments for mesh boundary flux.
!     "nc"    is number of spatial moments for mesh interior flux.
!     "nbd"   is a dimension of the coefficient matrices nb*nd

      PROGRAM SWEEP_ADAPTATIV

      USE SNQHRM
      USE SWEEP8ONE

      IMPLICIT NONE

      INTEGER, PARAMETER :: ng = 1, ndir = 3, ndim = 3
      INTEGER, PARAMETER :: nn = ndir*ng
      INTEGER, PARAMETER :: nc=4, nb=3, nd = ndir*8
      INTEGER, PARAMETER :: nbd = nb*ndim
      INTEGER, PARAMETER :: nx = 2, ny = 2, nz = 2
      INTEGER, PARAMETER :: nbfx = ny*nz,nbfy=nx*nz,nbfz=nx*ny
      INTEGER, PARAMETER :: nani=0,nhrm=1 
      INTEGER, PARAMETER :: nr = nx*ny*nz , nh = 1
      REAL, PARAMETER :: delt = 1.0
      INTEGER, PARAMETER :: nmat = 2
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

      INTEGER :: r,o
      
      zreg = 1
      DO r=5,nr
        zreg(r) = 2
      ENDDO

      sigs = 0.0
      sigt(:,1) = 0.15
      sigt(:,2) = 0.47
      
      bflx = 0.0
      bfly = 0.0
      bflz = 0.0

      bflx(:,1,:,:,:) = 1.2 
      bfly(:,1,:,:,:) = 1.2    
      bflz(:,1,:,:,:) = 1.2

      aflx = 0.0
      asrc = 0.0
      srcm = 0.0

      DO r=1,4
        srcm(:,r,1,1) = 0.18
      ENDDO
      DO r=5,nr
        srcm(:,r,1,1) = 0.564 
      ENDDO

      DO r = 1,nr
        DO o = 1,noct
            aflx(:,r,1,o) = 1.2
        ENDDO
      ENDDO


      print *,"av flux", aflx(:,:,1,1)
      
      
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
     &                     flxm)
        

      print *,"Flux final", aflx(:,:,1,1)
    !   print *, bflx  
    !   print *, bfly  
    !   print *, bflz  

      print *, 'FIN'
    !  
      END PROGRAM