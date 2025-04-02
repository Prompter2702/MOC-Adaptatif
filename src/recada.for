       SUBROUTINE SWEE3D_ADAPTIVE(
       ! inputs: dimensions 
     &                          nn,ng,nr,nh,nc,nmat,
     &                          nb,nbd,nbfx,nbfy,nbfz,
     &                          nx,ny,nz,
     &                          nani,nhrm,nd,ndir,
       ! inputs: cross section  
     &                          sigt,sigs,
       ! inputs: source and flux moments at previous iteration 
     &                          flxp,srcm,
       ! input/output: boundary flux 
     &                          bflx,bfly,bflz,
       ! input: angular quadrature & spherical harmonics
     &                          mu,eta,ksi,w,pisn,sphr,
       ! input: sing matrices to adapt coefficients to octants 
     &                          sgnc,sgni,sgne,sgnt,
       ! input: auxiliar memory for boundary angular fluxes, angular 
       ! source moments and self-scattering xs 
     &                          finc,fout,tmom,sigg,
       ! input: auxiliar memory angular flux & source 
     &                          aflx,asrc,dsrc,
       ! input: auxiliar memory for the interface angular flux 
     &                          xaux,yaux,zaux,
       ! input: auxiliar memory to drive angular mirror-reflection 
       ! or rotation/translation b.c.
     &                          rdir,
       ! input: pixel-to-medium array 
     &                          zreg, 
       ! input: pixel-to-medium array 
     &                          dira,dirf,
       ! input: logical to add angular source in case of
       ! time-dependent calculation 
     &                          lgki,
       ! output: flux moments 
     &                          flxm,
       ! Delta on each direction of the region
     &                          delt3)

      USE FLGOCT
      USE FLGBCD
      USE FLGINC
      USE SWEEP8ONE
      
      IMPLICIT NONE

      INTEGER, PARAMETER :: noct = 8, n8=8, ns=4
      
      INTEGER :: nn,ng,nr,nh,nc,nx,ny,nz,nmat
      INTEGER :: nb,nbd,nbfx,nbfy,nbfz
      INTEGER :: nani,nhrm,nd,ndir
      REAL    :: flxm(ng,nr,nh,nc)
      REAL    :: bflx(nn,nb,nbfx,2,noct),bfly(nn,nb,nbfy,2,noct),
     &           bflz(nn,nb,nbfz,2,noct)
      REAL    :: sigt(ng,nmat),sigs(ng,0:nani,nmat)
      REAL(KIND=8) :: mu(ndir),eta(ndir),ksi(ndir),w(ndir),pisn
      REAL    :: sphr(nhrm,nd)
      INTEGER :: sgnc(nc,nc,noct),sgni(nc,nbd,noct)
      INTEGER :: sgne(nbd,nc,noct),sgnt(nbd,nbd,noct)
      INTEGER :: zreg(nr)
      REAL    :: aflx(nn,nr,nc,noct),asrc(nn,nr,nc,noct)
      REAL    :: dsrc(nn,nr,nc,*)
      REAL    :: finc(nn,nbd),fout(nn,nbd)
      REAL    :: flxp(ng,nr,nh,nc),srcm(ng,nr,nh,nc),tmom(ng,nr,nh,nc)
      REAL    :: sigg(ng,nr)
      REAL    :: xaux(nn,nb),yaux(nn,nx,nb),zaux(nn,nx,ny,nb)
      INTEGER :: rdir(nd,*),dira(*),dirf(*)
      LOGICAL :: lgki
      REAL, INTENT(IN)     :: delt3(3)

      INTEGER :: xout,yout,zout,fst
      INTEGER :: oc,oct,d,dd,da,c,h,i

      REAL  :: aflx0(nn,nc), aflx1(nn,nc,n8)
      REAL  :: xshom0(ng), xshom1(ng,n8)
      REAL  :: asrcm0(ng,ndir,nc), asrcm1(ng,ndir,nc,n8)
      REAL  :: finc1(nn,nb,3,ns), fout1(nn,nb,3,ns)
      REAL  :: finc0(nn,nb,3), fout0(nn,nb,3)

      REAL  :: ccof(nn,nc,nc),icof(nn,nc,nbd)
      REAL  :: ecof(nn,nbd,nc),tcof(nn,nbd,nbd)
      REAL(KIND=8)  :: errcor(ng, ndir),errmul(ng, ndir)
      REAL     :: ccof8(nn,nc,nc,n8),icof8(nn,nc,nbd,n8)
      REAL     :: ecof8(nn,nbd,nc,n8),tcof8(nn,nbd,nbd,n8)
      LOGICAL  :: ok 
      INTEGER  :: ii,jj,kk,rr 

     

      CALL GMOM3D(ng,nani,nh,nc,nr,sigs,flxp,srcm,zreg,sigg,tmom)

!     Loop over octants of angular space in order defined by
!     octant list "olst3D".

      DO oc=8,8
         oct=olst3D(oc)
!        Index of outgoing side.
         xout=3-xinc(oct)
         yout=3-yinc(oct)
         zout=3-zinc(oct)
         
!        Sweep for all directions in octant "oct".
!        initial direction of the octant 

         da = (oct-1)*ndir+1

!        Directional source.

         CALL GIRSRC(nhrm,ng,ndir,nr,nh,nc,sphr(1,da),
     &               tmom(1,1,1,1),asrc(1,1,1,oct))
        
         IF(lgki)CALL SPXPY(nr*nc,dsrc(1,1,1,oct),asrc(1,1,1,oct))


!        Sweep.
        ! The boundary entries are computed on the exiting part
        ! in order to work only on the exiting part of the bdfl
        DO i=1,nbfx
            bflx(:,:,i,xout,oct) = bflx(:,:,i,xinc(oct),oct)
        ENDDO
        DO i=1,nbfy
            bfly(:,:,i,yout,oct) = bflx(:,:,i,yinc(oct),oct)
        ENDDO
        DO i=1,nbfz
            bflz(:,:,i,zout,oct) = bflx(:,:,i,zinc(oct),oct)
        ENDDO

!     Niveau 0 doit être déjà fait :
        CALL XSSRCHOMO0(nn,ng,nr,nx,ny, ndir,
     &                 1,nx,1,ny,1,nz,
     &                 sigt, zreg, aflx(:,:,:,oct), w,
     &                 asrc(:,:,:,oct), xshom0,
     &                 asrcm0)

        CALL ONE_COEF3D(nn,ndir,ng,mu,eta,ksi,
     &                delt3,xshom0,
     &                sgnc(:,:,oct),sgni(:,:,oct),
     &                sgne(:,:,oct),sgnt(:,:,oct),
     &                ccof,icof,ecof,tcof)
      

        CALL MERGEBOUND0(nn, ng,nb,
     &                   1,nx,1,ny,1,nz,
     &                   nx,ny,nz,
     &                   bflx(:,:,:,xinc,oct),
     &                   bfly(:,:,:,yinc,oct),
     &                   bflz(:,:,:,zinc,oct),
     &                   finc0)



        CALL SWEEP_ONEREGION(nn,2,asrcm0, finc0, aflx0, fout0,
     &                          ccof,icof,ecof,tcof)


        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc(oct), yinc(oct), zinc(oct),oct,
     &                             1,nx,1,ny,1,nz,
     &                             0,sigt, 
     &                             mu,eta,ksi,w,
     &                             sgnc(:,:,oct),sgni(:,:,oct),
     &                             sgne(:,:,oct),sgnt(:,:,oct),
     &                             zreg, asrc(:,:,:,oct),
     &                             aflx(:,:,:,oct),
     &                             bflx(:,:,:,xout,oct),
     &                             bfly(:,:,:,yout,oct),
     &                             bflz(:,:,:,zout,oct),
     &                             aflx0, aflx1,xshom0, xshom1,
     &                             asrcm0, asrcm1,finc0,finc1,
     &                             fout0,fout1,ccof,icof,ecof,tcof,
     &                             ccof8,icof8,ecof8,tcof8,
     &                             errcor,errmul,ok,delt3,ii,jj,kk,rr)  

         ! CALL SWEEP
!        Angular moments.
     
         DO c=1,nc
         DO h=1,nh
         fst = 0
         DO d = 1,ndir
            flxm(:ng,:nr,h,c) = flxm(:ng,:nr,h,c) + 
     &        (sphr(h,da)*w(d))*aflx(fst+1:fst+ng,:nr,c,oct)
            fst = fst + ng
         ENDDO
         ENDDO
         ENDDO
         
      ENDDO
     

!     Boundary conditions: specular reflection, albedo
    
      DO oc=1,8
         oct=olst3D(oc)

!        Index of outgoing side.

         xout=3-xinc(oct)
         yout=3-yinc(oct)
         zout=3-zinc(oct)
         
         CALL BCD2R3(oct,xout,yout,zout,
     &               bflx,bfly,bflz,mu,eta,
     &               ksi,w,pisn,rdir)
     
      ENDDO 
      
      END

!----------------------------------------------------------------------
!----------------------------------------------------------------------
      
      RECURSIVE SUBROUTINE RECADA_ONE_OCTANT(nn,ng,ndir,nr,nh,
     &                             nx,ny,nz,
     &                             xinc,yinc,zinc, oct,
     &                             imin,imax,jmin,jmax,kmin,kmax,niv,
     &                             sigt, 
     &                             mu,eta,ksi,w,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflx0, aflx1, xshom0, xshom1,
     &                             asrcm0, asrcm1,finc0,finc1,
     &                             fout0,fout1,ccof,icof,ecof,tcof,
     &                             ccof8,icof8,ecof8,tcof8,
     &                             errcor,errmul,ok,delt3,i,j,k,r)
     
      USE SWEEP8ONE
      USE SRCCORR 

      IMPLICIT NONE
     
      INTEGER, PARAMETER :: nc = 4, nb = 3, nbd = 9, n8 = 8,
     &                    ns = 4

      INTEGER, INTENT(IN) :: nn,ng,ndir,nr,nh,
     &           imin,imax,jmin,jmax,kmin,kmax,
     &           nx,ny,nz,niv,oct
      INTEGER, INTENT(IN) :: xinc,yinc,zinc
    !   REAL, INTENT(IN)    :: flxm(ng,nr,nh,nc)
      REAL, INTENT(IN)    :: sigt(ng,*)
      REAL(KIND=8), INTENT(IN) :: mu(ndir),eta(ndir),ksi(ndir),w(ndir)
      INTEGER, INTENT(IN) :: sgnc(nc,nc),sgni(nc,nbd)
      INTEGER, INTENT(IN) :: sgne(nbd,nc),sgnt(nbd,nbd)
      INTEGER, INTENT(IN) :: zreg(nr)
      REAL, INTENT(INOUT) :: asrc(nn,nr,nc)
      REAL, INTENT(INOUT) :: aflx(nn,nr,nc)
      REAL, INTENT(INOUT) :: bflx(nn,nb,ny*nz),bfly(nn,nb,nx*nz),
     &           bflz(nn,nb,nx*ny)
     
      REAL, INTENT(INOUT)  :: aflx0(nn,nc), aflx1(nn,nc,n8)
      REAL, INTENT(INOUT)  :: xshom0(ng), xshom1(ng,n8)
      REAL, INTENT(INOUT)  :: asrcm0(ng,ndir,nc), asrcm1(ng,ndir,nc,n8)
      REAL, INTENT(INOUT)  :: finc1(nn,nb,3,ns), fout1(nn,nb,3,ns)
      REAL, INTENT(INOUT)  :: finc0(nn,nb,3), fout0(nn,nb,3)

      REAL, INTENT(INOUT)  :: ccof(nn,nc,nc),icof(nn,nc,nbd)
      REAL, INTENT(INOUT)  :: ecof(nn,nbd,nc),tcof(nn,nbd,nbd)
      REAL(KIND=8),INTENT(INOUT)  :: errcor(ng, ndir),errmul(ng, ndir)
      REAL, INTENT(IN)     :: ccof8(nn,nc,nc,n8),icof8(nn,nc,nbd,n8)
      REAL, INTENT(IN)     :: ecof8(nn,nbd,nc,n8),tcof8(nn,nbd,nbd,n8)
      LOGICAL, INTENT(INOUT)  :: ok 
      REAL, INTENT(IN)        :: delt3(3)
      INTEGER, INTENT(INOUT)  :: i,j,k,r 


      REAL :: aflxtps1(nn,nc), asrcmtps1(ng,ndir,nc), xshomtps1(ng)
      REAL :: aflxtps2(nn,nc), asrcmtps2(ng,ndir,nc), xshomtps2(ng)
      REAL :: aflxtps3(nn,nc), asrcmtps3(ng,ndir,nc), xshomtps3(ng)
      REAL :: aflxtps4(nn,nc), asrcmtps4(ng,ndir,nc), xshomtps4(ng)
      REAL :: aflxtps5(nn,nc), asrcmtps5(ng,ndir,nc), xshomtps5(ng)
      REAL :: aflxtps6(nn,nc), asrcmtps6(ng,ndir,nc), xshomtps6(ng)
      REAL :: aflxtps7(nn,nc), asrcmtps7(ng,ndir,nc), xshomtps7(ng)
      REAL :: aflxtps8(nn,nc), asrcmtps8(ng,ndir,nc), xshomtps8(ng)
      
      fout1 = 0.0
! Error check by source-correction estimation

      CALL SRCCOR(ng,ndir,delt3/(2**niv),mu,eta,ksi,asrcm0,
     &            xshom0,errcor)

!     7/ query ok/ko? 
!     Définir un critère 
      IF(imax-imin>0 .AND. jmax-jmin>0 .AND. kmax-kmin>0) THEN
        ok = .FALSE.
      ELSE
        ok = .TRUE.
      ENDIF

      IF (.NOT. ok) THEN 

        CALL XSSRCHOMO1(nn,ng,nr,nx,ny, ndir,
     &                 imin, imax,jmin,jmax,kmin,kmax,
     &                 sigt, zreg, aflx, w, asrc, xshom1,
     &                 asrcm1)

        CALL MERGEBOUND1(nn, ng,nb,
     &                   imin, imax,jmin,jmax,kmin,kmax,
     &                   nx,ny,nz,
     &                   bflx, bfly, bflz, finc1)


        CALL EIGHT_COEF3D(nn,ndir,ng,mu,eta,ksi,
     &                    delt3/(2**(niv+1)),xshom1,
     &                    sgnc,sgni,sgne,sgnt,
     &                    ccof8,icof8,ecof8,tcof8)


        CALL SWEEP_8REGIONS(nn,2,asrcm1, finc1, aflx1, fout1,
     &                     ccof8,icof8,ecof8,tcof8, xinc, yinc, zinc)


    !   CALL SRC2LVL(ng, ndir, srcm1, sigt, errmul)
      ENDIF

      IF (ok) THEN 
        ! - write the angular flux of node-0 on the fine-mesh flux 
        DO i=imin,imax
            DO j=jmin,jmax
                DO k=kmin,kmax
                    r = ((k-1)*ny + (j-1))*nx + i
                    aflx(:,r,:) = aflx0(:,:)
                ENDDO
            ENDDO
        ENDDO
!      - write the outgoing angular flux of node-0 
!            on the fine-mesh boundary flux
         CALL SPLITBOUND0(nn, ng,nb,
     &                   imin, imax,jmin,jmax,kmin,kmax,
     &                   nx,ny,nz,
     &                   bflx, bfly, bflz, fout0)
        
          !         - return 
          RETURN

      ELSE IF((imax-imin)>1) THEN
!         15/ si ko : 
!        - call recursively the subroutine for the 8 nodes of level-1


        asrcmtps1 = asrcm1(:,:,:,xinc+2*(yinc-1) + 4*(zinc-1))  
        aflxtps1 = aflx1(:,:, xinc+2*(yinc-1) + 4*(zinc-1))
        xshomtps1 = xshom1(:,  xinc+2*(yinc-1) + 4*(zinc-1))

        asrcmtps2 = asrcm1(:,:,:,3-xinc+2*(yinc-1) + 4*(zinc-1))  
        aflxtps2 = aflx1(:,:, 3-xinc+2*(yinc-1) + 4*(zinc-1))
        xshomtps2 = xshom1(:,  3-xinc+2*(yinc-1) + 4*(zinc-1))

        asrcmtps3= asrcm1(:,:,:,xinc+2*(2-yinc) + 4*(zinc-1))  
        aflxtps3= aflx1(:,:, xinc+2*(2-yinc) + 4*(zinc-1))
        xshomtps3= xshom1(:,  xinc+2*(2-yinc) + 4*(zinc-1))

        asrcmtps4= asrcm1(:,:,:,3-xinc + 2*(2-yinc) + 4*(zinc-1))  
        aflxtps4= aflx1(:,:, 3-xinc + 2*(2-yinc) + 4*(zinc-1))
        xshomtps4= xshom1(:,  3-xinc + 2*(2-yinc) + 4*(zinc-1))

        asrcmtps5= asrcm1(:,:,:,xinc + 2*(yinc-1) + 4*(2-zinc))  
        aflxtps5= aflx1(:,:, xinc + 2*(yinc-1) + 4*(2-zinc))
        xshomtps5= xshom1(:,  xinc + 2*(yinc-1) + 4*(2-zinc))

        asrcmtps6= asrcm1(:,:,:,3-xinc + 2*(yinc-1) + 4*(2-zinc))  
        aflxtps6= aflx1(:,:, 3-xinc + 2*(yinc-1) + 4*(2-zinc))
        xshomtps6= xshom1(:,  3-xinc + 2*(yinc-1) + 4*(2-zinc))  

        asrcmtps7= asrcm1(:,:,:,xinc+2*(2-yinc) + 4*(2-zinc))  
        aflxtps7= aflx1(:,:, xinc+2*(2-yinc) + 4*(2-zinc))
        xshomtps7= xshom1(:,  xinc+2*(2-yinc) + 4*(2-zinc))

        asrcmtps8= asrcm1(:,:,:,3-xinc + 2*(2-yinc) + 4*(2-zinc))  
        aflxtps8= aflx1(:,:, 3-xinc + 2*(2-yinc) + 4*(2-zinc))
        xshomtps8= xshom1(:,  3-xinc + 2*(2-yinc) + 4*(2-zinc))

        asrcm0 = asrcmtps1
        aflx0  = aflxtps1
        xshom0 = xshomtps1

        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &   (2*imin + (xinc-1)*(imax-imin+2))/2 ,
     &   (imin + imax + (xinc-1)*(imax-imin))/2,
     &   (2*jmin + (yinc-1)*(jmax-jmin+2))/2,
     &   (jmin + jmax + (yinc-1)*(jmax-jmin))/2,
     &   (2*kmin + (zinc-1)*(kmax-kmin+2))/2,
     &   (kmin + kmax + (zinc-1)*(jmax-jmin))/2,
     &                             niv+1,
     &                             sigt,
     &                             mu,eta,ksi,w,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflx0, aflx1,xshom0, xshom1,
     &                             asrcm0, asrcm1,finc0,finc1,
     &                             fout0,fout1,ccof,icof,ecof,tcof,
     &                             ccof8,icof8,ecof8,tcof8,
     &                             errcor,errmul,ok,delt3,i,j,k,r)
      
      asrcm0 = asrcmtps2
      aflx0  = aflxtps2
      xshom0 = xshomtps2

        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &   (2*imin + (2-xinc)*(imax-imin+2))/2,
     &   (imin + imax + (2-xinc)*(imax-imin))/2,
     &   (2*jmin + (yinc-1)*(jmax-jmin+2))/2,
     &   (jmin + jmax + (yinc-1)*(jmax-jmin))/2,
     &   (2*kmin + (zinc-1)*(kmax-kmin+2))/2,
     &   (kmin + kmax + (zinc-1)*(jmax-jmin))/2,
     &                             niv+1,
     &                             sigt,
     &                             mu,eta,ksi,w,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflx0, aflx1,xshom0, xshom1,
     &                             asrcm0, asrcm1,finc0,finc1,
     &                             fout0,fout1,ccof,icof,ecof,tcof,
     &                             ccof8,icof8,ecof8,tcof8,
     &                             errcor,errmul,ok,delt3,i,j,k,r)

      asrcm0 = asrcmtps3
      aflx0  = aflxtps3
      xshom0 = xshomtps3

        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &   (2*imin + (xinc-1)*(imax-imin+2))/2 ,
     &   (imin + imax + (xinc-1)*(imax-imin))/2,
     &   (2*jmin + (2-yinc)*(jmax-jmin+2))/2,
     &   (jmin + jmax + (2-yinc)*(jmax-jmin))/2,
     &   (2*kmin + (zinc-1)*(kmax-kmin+2))/2,
     &   (kmin + kmax + (zinc-1)*(jmax-jmin))/2,
     &                             niv+1,
     &                             sigt,
     &                             mu,eta,ksi,w,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflx0, aflx1,xshom0, xshom1,
     &                             asrcm0, asrcm1,finc0,finc1,
     &                             fout0,fout1,ccof,icof,ecof,tcof,
     &                             ccof8,icof8,ecof8,tcof8,
     &                             errcor,errmul,ok,delt3,i,j,k,r)

      asrcm0 = asrcmtps4
      aflx0  = aflxtps4
      xshom0 = xshomtps4

        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &   (2*imin + (2-xinc)*(imax-imin+2))/2,
     &   (imin + imax + (2-xinc)*(imax-imin))/2,
     &   (2*jmin + (2-yinc)*(jmax-jmin+2))/2,
     &   (jmin + jmax + (2-yinc)*(jmax-jmin))/2,
     &   (2*kmin + (zinc-1)*(kmax-kmin+2))/2,
     &   (kmin + kmax + (zinc-1)*(jmax-jmin))/2,
     &                             niv+1,
     &                             sigt,
     &                             mu,eta,ksi,w,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflx0, aflx1,xshom0, xshom1,
     &                             asrcm0, asrcm1,finc0,finc1,
     &                             fout0,fout1,ccof,icof,ecof,tcof,
     &                             ccof8,icof8,ecof8,tcof8,
     &                             errcor,errmul,ok,delt3,i,j,k,r)


!-----------------------------------------------------------------------
      asrcm0 = asrcmtps5
      aflx0  = aflxtps5
      xshom0 = xshomtps5
        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &   (2*imin + (xinc-1)*(imax-imin+2))/2 ,
     &   (imin + imax + (xinc-1)*(imax-imin))/2,
     &   (2*jmin + (yinc-1)*(jmax-jmin+2))/2,
     &   (jmin + jmax + (yinc-1)*(jmax-jmin))/2,
     &   (2*kmin + (2-zinc)*(kmax-kmin+2))/2,
     &   (kmin + kmax + (2-zinc)*(kmax-kmin))/2,
     &                             niv+1,
     &                             sigt,
     &                             mu,eta,ksi,w,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflx0, aflx1,xshom0, xshom1,
     &                             asrcm0, asrcm1,finc0,finc1,
     &                             fout0,fout1,ccof,icof,ecof,tcof,
     &                             ccof8,icof8,ecof8,tcof8,
     &                             errcor,errmul,ok,delt3,i,j,k,r)

      asrcm0 = asrcmtps6
      aflx0  = aflxtps6
      xshom0 = xshomtps6


        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &   (2*imin + (2-xinc)*(imax-imin+2))/2,
     &   (imin + imax + (2-xinc)*(imax-imin))/2,
     &   (2*jmin + (yinc-1)*(jmax-jmin+2))/2,
     &   (jmin + jmax + (yinc-1)*(jmax-jmin))/2,
     &   (2*kmin + (2-zinc)*(kmax-kmin+2))/2,
     &   (kmin + kmax + (2-zinc)*(kmax-kmin))/2,
     &                             niv+1,
     &                             sigt,
     &                             mu,eta,ksi,w,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflx0, aflx1,xshom0, xshom1,
     &                             asrcm0, asrcm1,finc0,finc1,
     &                             fout0,fout1,ccof,icof,ecof,tcof,
     &                             ccof8,icof8,ecof8,tcof8,
     &                             errcor,errmul,ok,delt3,i,j,k,r)      
      asrcm0 = asrcmtps7
      aflx0  = aflxtps7
      xshom0 = xshomtps7
        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &   (2*imin + (xinc-1)*(imax-imin+2))/2 ,
     &   (imin + imax + (xinc-1)*(imax-imin))/2,
     &   (2*jmin + (2-yinc)*(jmax-jmin+2))/2,
     &   (jmin + jmax + (2-yinc)*(jmax-jmin))/2,
     &   (2*kmin + (2-zinc)*(kmax-kmin+2))/2,
     &   (kmin + kmax + (2-zinc)*(kmax-kmin))/2,
     &                             niv+1,
     &                             sigt,
     &                             mu,eta,ksi,w,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflx0, aflx1,xshom0, xshom1,
     &                             asrcm0, asrcm1,finc0,finc1,
     &                             fout0,fout1,ccof,icof,ecof,tcof,
     &                             ccof8,icof8,ecof8,tcof8,
     &                             errcor,errmul,ok,delt3,i,j,k,r)       
      asrcm0 = asrcmtps8
      aflx0  = aflxtps8
      xshom0 = xshomtps8
        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &   (2*imin + (2-xinc)*(imax-imin+2))/2,
     &   (imin + imax + (2-xinc)*(imax-imin))/2,
     &   (2*jmin + (2-yinc)*(jmax-jmin+2))/2,
     &   (jmin + jmax + (2-yinc)*(jmax-jmin))/2,
     &   (2*kmin + (2-zinc)*(kmax-kmin+2))/2,
     &   (kmin + kmax + (2-zinc)*(kmax-kmin))/2,
     &                             niv+1,
     &                             sigt,
     &                             mu,eta,ksi,w,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflx0, aflx1,xshom0, xshom1,
     &                             asrcm0, asrcm1,finc0,finc1,
     &                             fout0,fout1,ccof,icof,ecof,tcof,
     &                             ccof8,icof8,ecof8,tcof8,
     &                             errcor,errmul,ok,delt3,i,j,k,r)       

    !- write the outgoing angular flux of node-0 on the fine-mesh
    ! boundary flux
         RETURN
      
        ELSE 
!   Case where the error is not guaranted but the pixel-lvl is reached
            CALL SPLITAFLX1(nn,nr,nc,nx,ny,nz,
     &                      imin,imax,jmin,jmax,kmin,kmax,
     &                      aflx,aflx1)


            CALL SPLITBOUND1(nn,ng,nb,
     &                      imin, imax,jmin,jmax,kmin,kmax,
     &                      nx,ny,nz,
     &                      bflx, bfly, bflz, fout1)

      ENDIF

      END SUBROUTINE RECADA_ONE_OCTANT