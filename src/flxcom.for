        MODULE FLXCOM

        ! Dimensions of arrays and flags used in flux solver.

        ! Dimensions

        INTEGER :: ndim,noct      ! Number of dimensions and octants
                                ! (geometry dependent parameters).

        INTEGER :: ng,nfg         ! Number of groups, fast groups.
        INTEGER :: nm             ! Number of media.
        INTEGER :: na             ! Max anisotropy order.
        INTEGER :: nani           ! Max anisotropy order for flux solver.
        INTEGER :: nzon           ! Number of zones.
        INTEGER :: nasg           ! Number of Zone to medium assignments.
        INTEGER :: nzx,nzy,nzz    ! Dimensions of zone mesh.
        INTEGER :: nx,ny,nz       ! Dimensions of computational mesh.
        INTEGER :: npx,npy,npz    ! (/nx+1,ny+1,nz+1/)
        INTEGER :: nr,np          ! Number of regions, points.
        INTEGER :: nbfx,nbfy,nbfz ! Number of meshes on outer surfaces.
        INTEGER :: nbix,nbiy,nbiz ! Number of meshes interfaces.
        INTEGER :: nbf(3)         ! (/nbfx,nbfy,nbfz/)
        INTEGER :: nbi(3)         ! (/nbix,nbiy,nbiz/)

        LOGICAL :: gdia           ! .T. if diagonal symmetry (unfolded)

        INTEGER :: nfr            ! Sum over regions of number of
                                ! fissile isotopes.
        INTEGER :: ncgs           ! Number of coarse group for fission source & spectrum

        INTEGER :: nd,ndir        ! Total number of discrete directions
                                ! and number of directions per octant.

        INTEGER :: nh             ! Number of flux harmonics.
        INTEGER :: nhrm           ! Number of spherical harmonics coeff.
        INTEGER :: nb,nc          ! Number of spatial flux moments on
                                ! boundary interfaces and in interior

        INTEGER :: nbd            ! = nb * ndim
        INTEGER :: nrid           ! Number of mesh cell coefficients.
        !38

        ! Boundary conditions

        INTEGER :: typc(2,3)      ! b.c. flags.
        LOGICAL :: lgsp(2,3)      ! .T. if specular reflection.
        LOGICAL :: lgtr(3)        ! .T. if translation.
        LOGICAL :: lgal(2,3)      ! .T. if albedo condition.

        LOGICAL :: lgfr           ! .T. if full reflection without
                                ! transfer on all boundaries.

        LOGICAL :: lgir           ! .T. if specular reflection on all
                                ! starting boundaries.

        LOGICAL :: lgbs(2,3)      ! .T. if boundary source
        !7

        ! Type of problem

        LOGICAL :: lgei           ! Eigenvalue problem if .TRUE.
        LOGICAL :: lgki           ! Kinetic problem if .TRUE.
        LOGICAL :: lgno           ! noise calculation if .TRUE.
        LOGICAL :: lgfs           ! include fissions if .TRUE.

        ! Calculational options

        INTEGER :: meth           ! Method of solution, 0=diamond,
                                ! 1=nodal, 2=characteristics.

        INTEGER :: mord           ! Flux expansion flag,
                                ! 0=diamond, 1=constant, 2=linear,
                                ! 3=bilinear, 4=parabolic.

        INTEGER :: no             ! Nodal expansion flag if NEM,
                                ! 0=parabolic, 2=quartic.

        INTEGER :: trcr           ! Transport correction flag
                                ! 0=none, 1=diagonal, 2=Bell-Hansen,
                                ! 3=Cesaro

        ! Iterations, convergence

        INTEGER :: nfgs,ntgs      ! Number of groups treated as fast,
                                ! number of thermal groups to be solved
                                ! simultaneously.

        INTEGER :: nout           ! Max. number of outer iterations.
        INTEGER :: ninn,nin1,nin2 ! Max. number of inner iterations.
        INTEGER :: nthr,nth1,nth2 ! Max. number of thermal iterations.
        INTEGER :: nacc           ! Max. number of acc. eq. iterations.

        INTEGER :: icco  ! Error norm flag for the convergence check


        REAL    :: epso,epsf,epsi ! Precision on eigenvalue, fission
        REAL    :: epst,epsa,epsb ! integral and inners, thermal and
                                ! acceleration iterations. and boundary flx

        INTEGER :: acco           ! Type of acceleration for outers
                                ! 0=none, 1=DSA, 2=Chebyshev.

        INTEGER :: acct           ! Type of acceleration for thermals
                                ! 0=none, 1=group rescaling.

        INTEGER :: acci           ! Type of acceleration for inners
                                ! 0=none, 1=DSA, 2=BPA.

        INTEGER :: acis           ! Type of discretization scheme for
                                ! acc. equation of inners
                                ! 0=diamond consistent, 1=wdd.

        INTEGER :: acim           ! Metod of solution for
                                ! acc. equation of inners,
                                ! 1=alternating directions, 2=multigrid.

        INTEGER :: mgri           ! Max. number of grids for multigrid
                                ! solver.

        INTEGER :: ngrd           ! Number of grids for multigrid.
        INTEGER :: nhmg           ! Length of multigrid mesh array.

        INTEGER :: nche           ! Number of free iterations before
                                ! Chebyshev acceleration.

        LOGICAL :: lgac           ! Acceleration of inners if .T.
        INTEGER :: inif           ! Unit of initial flux guess file.

        LOGICAL :: lgka           ! Angular flux kept if .T.
        LOGICAL :: lgkc           ! Interface currents kept if .T.
        INTEGER :: inda           ! Increment of direction count for
                                ! angular flux storage (0 or 1).
        INTEGER :: indc           ! Case flag for interface current
                                ! storage (0 or 1).

        INTEGER      :: tinn,iout ! Total number of inners, outers.
        INTEGER      :: nine      ! Number of inners per outer.
        INTEGER      :: tthe      ! Total number of thermals.
        REAL         :: taci      ! Total number of acc. eq. iterations.
        REAL(KIND=8) :: timi,timo ! CPU time in inners, outers.
        REAL(KIND=8) :: timc      ! CPU time for coefficients (Transport+Acceleration) .
        REAL         :: tmem      ! Total auxiliary memory.
        !48

        ! Acceleration

        INTEGER :: nacm           ! Max flux harmonic number to be
                                ! accelerated.

        INTEGER :: ndfc,ndsa      ! Number of anisotropy components for
                                ! DSA stuff (depends on number of
                                ! moments to be accelerated).

        INTEGER :: nwde,nwds      ! Number of coefficients for WDD acc.
        INTEGER :: nwdg,nwdc      ! schemes.

        ! Printout flags

        LOGICAL :: lgeh           ! Echo on console if .T.
        LOGICAL :: lgou,lgin,lgth ! Iteration status for outers, inners.

        INTEGER :: ctrp           ! Common transfer profile if 1.
        REAL    :: cmem           ! Memory for coefficients in words.
        INTEGER :: gstc           ! Number of groups for coefficients.
        !13


        ! Adjoint calculation flag

        LOGICAL :: adjo

        ! Auxiliary index for outer iteration
        !    (auxiliary stuff must be elsewhere)

        INTEGER :: pass

        ! parameter and index to perform multigroup sweeping for both
        ! direct and adjoint calculations.

        INTEGER  :: nfgfr,nfgls,pasfg ! first and last fast group and jump
        INTEGER  :: nthfr,nthls,pasth ! first and last fast group and jump

        INTEGER , DIMENSION(3) :: itrlst ! iteration block list
        
        ! Multigroup Synthetic Acceleration parameters :
        INTEGER :: ntia ! n. of thermals iteration for Synthetic Thermal Iterations
        INTEGER :: noia ! n. of outers iteration for Synthetic Outer Iterations
        REAL    :: epta ! precision for Synthetic Thermal Iterations (error mesured on flux integral)
        REAL    :: epoa ! precision for Synthetic Outer Iterations (error mesured on flux integral)
        !13

        ! BPA stuff:
        ! ---------
        INTEGER :: naacc
        INTEGER :: nhacc
        INTEGER :: nhha

        ! BPA for HCC
        INTEGER :: nbsm !  n. of incident boundary surfaces (per octant)
        INTEGER :: ncbpa,nibpa,nebpa,ntbpa ! matrix dimensions for BPA coefficients

        ! HCC stuff
        ! ---------
        LOGICAL :: lghc
        LOGICAL :: lgcmem ! =.TRUE. if HCCs coefficients are stored in memory
        INTEGER :: npin ! total n. of HCC
        INTEGER :: nhh ! dimension of simmetric collision matrix w.r.t. the angular moments, nhh = nh*(nh+1)/2
        ! max spatial dimensions dimension for a HCC
        INTEGER :: nrmax ! max n. of regions in a HCC
        INTEGER :: nimax,nomax ! max n. of incoming/outgoing surfaces (over all angles) in a HCC
        INTEGER :: nsmax ! max n. of surfaces in a HCC
        ! sweeping data
        INTEGER :: nfrt  ! max n. of fronts (for all directions)
        INTEGER :: ntrk  ! n. of tracking types (geometrical types)
        ! geometry data
        INTEGER :: nbox  ! n. of boxes
        INTEGER :: nreg  ! total n. of regions in the geometry
        INTEGER :: nsur ! total n. of surfaces in the geometry
        ! coefficients dimensions
        INTEGER :: ncc,nch,ncb,nbb ! nch=nc*nc*nhh ; ncb=nc*nb ; nbb=nb*nb
        INTEGER :: ndofc,ndofi,ndofe,ndoft ! max n. of volume/incoming/outgoing/surface d.o.f.
        INTEGER :: cdim,idim,edim,tdim ! total number of elements for collision/incoming/escape/transfer matrix
        INTEGER :: nsub(3),nsxa(3) ! n. of subdivision per axis & n. of surfaces per axis
        INTEGER :: nsubx,nsxax ! max ... "
        INTEGER :: ncofe,ncofi,ncoft ! n. of sign-matrix types
        INTEGER :: nunt ! n. of unit-files for coefficients (=0 then all in memory)
        INTEGER :: mxpk ! max n. of groups
        INTEGER :: ngpk ! max n. of groups for coefficients stored in a file 
        LOGICAL :: lghco ! = .TURE. if the collision matrix is stored per angular moment (integrated over the angle)
                        ! = .FALSE. is stored per angle
        LOGICAL :: lgprtb ! = .TRUE. if the peturbation strategy is set on for limintg coefficient computation
        
        ! data for TRIANGLEs
        
        INTEGER :: ncmax !maximum dimension of the unknown vector
        INTEGER :: ncmin !minimum dimension of the unknown vector
        !INTEGER :: ncof  !total number of coefficients in the coefficient matrices
        !cgfn provvisoriamente come parametro
        INTEGER , PARAMETER :: ncof = 35
        INTEGER :: nbnd  !total number of boundaries
        INTEGER :: nbinx,nboux ! max number of incoming/outgoing boundaries
        !43

        LOGICAL :: lgddm
        LOGICAL :: lgomp
        LOGICAL :: acco_RHCG

        ! 46 
        INTEGER , DIMENSION(2) :: fsi_cpp_ptr
        LOGICAL :: lgap3
        LOGICAL :: doAutop
        LOGICAL :: doAutopOC
        LOGICAL :: doCofOnLine

        ! Multi-Pn method variables
        LOGICAL :: mpnf           ! Logical for the MPN method (mpnf.EQ..T.)
        LOGICAL :: mpnnc          ! Logical for the MPN not consistent mode (.T. not consitent enablad)
                                !                                         (.F. consitent - default)
        INTEGER :: mpne           ! Solid-angle discretization option:
                                ! 1=Rectangular shaped element
                                ! 2=Triangular shaped element 
        INTEGER :: mpno           ! Local angular expansion order:
                                ! 1=Constant, 2=Linear
                                ! 3=Bilinear, 4=Quadratic
                                ! 5=Cubic
        INTEGER :: mpnsp          ! IF mpne=1: number of solid-angle elements
                                ! along the polar direction, per quadrant.
        INTEGER :: mpnsa          ! IF mpne=1: number of solid-angle elements
                                ! along the azimuthal direction, per quadrant.
                                ! IF mpne=2: number of nodes of the solid-angle
                                ! triangular discretization, per quadrant.
        INTEGER :: mpnqo          ! MPn quadrature rule order per solid angle.
        INTEGER :: mpnqt          ! MPN quadrature type:
                                !  IF (mpne=1):
                                !    - 1=UNIFORM
                                !    - 2=GaussLegendre
                                !    - 3=ChebyshevTimesLegendre
                                !  IF (mpne=2):
                                !    - 1=uniform (centers of a trianglular decomposition)
        INTEGER :: nnode
        INTEGER (KIND=4) :: nsa   ! Number of solid angles   
        INTEGER (KIND=4) :: nfq   ! Number of quadrature points per solid-angle

        INTEGER :: nhc            ! Number of angular base functions
                                ! of the flux surface expansion 
        INTEGER :: nc_nh !=nc*nh
        INTEGER :: nb_nhc !=nb*nhc
        INTEGER :: nbd_nhc !=nbd*nhc
        INTEGER :: ncm            ! Characteristic size for MPN matrixes construction
                                ! ncm=nc if not conformal, ncm=nc_nh if conformal

        !FCS logical, Local to the subdomain solver
        LOGICAL :: lgfcs          ! True when the FCS has been calculated
                                ! This adds FCS to the source (srcm) and
                                ! adds the uncollided flux to the flux
                                ! once the collided-flux calculation is over.

        
    !$OMP  THREADPRIVATE(  
    !$OMP1 ndim ,noct ,ng   ,nfg   ,nm    ,na    ,nani  ,nzon  ,nasg  ,nzx ,
    !$OMP1 nzy  ,nzz  ,nx   ,ny    ,nz    ,npx   ,npy   ,npz   ,nr    ,np  ,
    !$OMP1 nbfx ,nbfy ,nbfz ,nbix  ,nbiy  ,nbiz  ,nbf   ,nbi   ,nfr   ,nd  ,
    !$OMP1 ndir ,nh   ,nhrm ,nb    ,nc    ,nbd   ,nrid  ,typc  ,lgsp  ,lgtr,
    !$OMP1 lgal ,lgfr ,lgir ,lgbs  ,lgei  ,meth  ,mord  ,trcr  ,nfgs  ,ntgs,
    !$OMP1 nout ,ninn ,nin1 ,nin2  ,nthr  ,nth1  ,nth2  ,nacc  ,epso  ,epsf,
    !$OMP1 epsi ,epst ,epsa ,acco  ,acct  ,acci  ,acis  ,acim  ,mgri  ,ngrd,
    !$OMP1 nhmg ,nche ,lgac ,inif  ,lgka  ,lgkc  ,inda  ,indc  ,tinn  ,iout,
    !$OMP1 nine ,tthe ,taci ,timi  ,timo  ,timc  ,tmem  ,nacm  ,ndfc  ,ndsa,
    !$OMP1 nwde ,nwds ,nwdg ,nwdc  ,lgeh  ,lgou  ,lgin  ,lgth  ,ctrp  ,cmem,
    !$OMP1 gstc, adjo ,pass ,nfgfr,nfgls ,pasfg ,nthfr ,nthls ,pasth ,itrlst,
    !$OMP1 ntia, noia ,epta ,epoa ,naacc ,nhacc ,nhha  ,nbsm  ,lgfs,
    !$OMP1 ncbpa,nibpa,nebpa,ntbpa ,lghc  ,lgcmem,npin  ,nhh ,
    !$OMP1 nrmax,nimax,nomax,nsmax ,nfrt  ,ntrk  ,nbox  ,nreg  ,nsur  ,ncc ,
    !$OMP1 nch  ,ncb  ,nbb  ,ndofc ,ndofi ,ndofe ,ndoft ,cdim  ,idim  ,edim,
    !$OMP1 tdim ,nsub ,nsxa ,nsubx ,nsxax ,ncofe ,ncofi ,ncoft ,nunt  ,mxpk,
    !$OMP1 ngpk ,lghco,lgprtb,lgddm,ncgs  ,no    ,lgomp ,acco_RHCG,
    !$OMP1 fsi_cpp_ptr,lgap3,doAutop,doAutopOC,doCofOnLine,
    !$OMP1 icco , gdia , lgno, lgki, epsb,
    !$OMP1 mpnf, mpnnc, mpne, mpno, mpnsp, mpnsa, mpnqo, mpnqt, 
    !$OMP1 nnode, nsa, nfq, nhc, nc_nh, nb_nhc, nbd_nhc, ncm , lgfcs)

        END MODULE
