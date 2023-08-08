CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                          C
C          SUBROUTINE 8-NODED SOLID SHELL ELEMENT          C
C                                                          C   
C   Author: Pavan Kumar Asur , Jose Reinoso                C
C  pavan.asur@imtlucca.it                                  C
C            Jan 2022                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      module kvisual
         IMPLICIT NONE
      real*8 USRVAR(70000,11,8), Ac
      integer nelem
	  INTEGER,PARAMETER :: allelem = 4452
      save
      end module
!-----------------------------------------------------------------------
!------------------------------------------------------------------------                                    
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3 LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)

      use kvisual
      INCLUDE 'ABA_PARAM.INC'
         !IMPLICIT NONE

      ! ----------------------------------------------------------
      ! variables defined in uel, passed back to Abaqus
      ! ----------------------------------------------------------
       REAL*8 rhs(ndofel,*),amatrx(ndofel,ndofel),svars(*),energy(8),
     +  pnewdt
       REAL*8 U(NDOFEL),DU(NDOFEL,*),V(NDOFEL),A(NDOFEL)
 
      ! ----------------------------------------------------------
      ! variables passed into UEL
      ! ----------------------------------------------------------
      INTEGER ndofel,nrhs,nsvars,nprops,mcrd,nnode,jtype,kstep,kinc,
     +  jelem,ndload,jdltyp(mdload,*),npredf,lflags(*),mlvarx,mdload,
     +  jprops(*),njprop
      !
      REAL*8 props(*),coords(mcrd,nnode),time(2),dtime,params(*),
     +  adlmag(mdload,*),predef(2,npredf,nnode),ddlmag(mdload,*),period
      ! ----------------------------------------------------------

!      INTEGER,PARAMETER :: nnode = 8 !number of element nodes
      INTEGER,PARAMETER :: ndim  = 3 !dimension
      INTEGER,PARAMETER :: PFTol  = 0.5d0 !CHange this if you want to change the threshold for EAS and ANS
      INTEGER,PARAMETER :: ntens = 6 !number of stress components
      INTEGER,PARAMETER :: ndof  = 4 ! ndof per node
      
      INTEGER,PARAMETER :: ngauss = 2 !number of integration points(IP) in-plane
      INTEGER,PARAMETER :: ngausst = 2 !number of IP over thickness
      
      INTEGER,PARAMETER :: nEAS = 7 !number of EAS parameters
      
      INTEGER,PARAMETER :: nANSshear = 2 !number of ANS in shear
      INTEGER,PARAMETER :: nANStrap = 4 !number of ANS due to trapizoidal
      INTEGER,PARAMETER :: nsvarsip = 5
	        REAL*8 statevLocal(nsvarsip)
      INTEGER isvinc
      !----------------------------------------------------
      integer i,j,k,m,kinkt,ik,jk,ii     
      REAL*8 xref(ndim,nnode) !REFERENCE CONFIG COORDINATES
      REAL*8 xcur(ndim,nnode) !CURRENT CONFIG COORDINATES      
      real*8 DeltaU(ndim*nnode,1)  
      !Weights and position for integration
      !-----------------------------------
      REAL*8 posgp(ngauss),weightgp(ngauss) !Position and weight IP
      REAL*8 posgpt(ngausst),weightgpt(ngausst) !Position and weight IP in thickness direction
      REAL*8 eta,phi, weta, wphi, xi, wxi 
      !-----------------------------------
      REAL*8 shapef(nnode),dshape(ndim,nnode)
      !-----------------------------------
      REAL*8 gkovr(ndim,ndim),gkonr(ndim,ndim)
      REAL*8 gmkovr(ndim,ndim),gmkonr(ndim,ndim), detr   
          
      REAL*8 gkovc(ndim,ndim),gkonc(ndim,ndim)
      REAL*8 gmkovc(ndim,ndim),gmkonc(ndim,ndim), detc 
      !-----------------------------------
      REAL*8 Bop(ntens,nnode*ndim), Jacobtr
      !---------------------------------
      REAL*8 Cmat(ntens,ntens), stress(ntens), strain(ntens)
      !--------------------------------
      REAL*8 emod,enu,matflag,KFreeEner,kFreeEner_pos,kFreeEner_neg
      !--------------------------------
      REAL*8 Kdd(nnode*ndim,nnode*ndim) !stiffness matrix at IP level
      REAL*8 Fd(nnode*ndim)  !internal force at IP level
      REAL*8 Kgdd(nnode*ndim,nnode*ndim)
      !--------------------------------
      REAL*8 Kedd(nnode*ndim,nnode*ndim) !stiffness matrix at element level
      REAL*8 Fed(nnode*ndim) !internal force at element level
      REAL*8 Kegdd(nnode*ndim,nnode*ndim),NLGEOMflag 
!---------------------------------------------------------         
      !--------------------------------
      !--  EAS MAGNITUDES
      !--------------------------------
      REAL*8 EASpar
      REAL*8 shapef0(nnode)
      REAL*8 dshape0(ndim,nnode)
        REAL*8 gkovr0(ndim,ndim)
        REAL*8 gkonr0(ndim,ndim)
        REAL*8 gmkovr0(ndim,ndim)
        REAL*8 gmkonr0(ndim,ndim)
        REAL*8 detr0
        REAL*8 gkovc0(ndim,ndim)
        REAL*8 gkonc0(ndim,ndim)
        REAL*8 gmkovc0(ndim,ndim)
        REAL*8 gmkonc0(ndim,ndim)
        REAL*8 detc0
!-------------------------------- 
        REAL*8 Meas(6,nEAS) !30 is the maximum number of EAS paramaters allowed
        REAL*8 Weasloc(6,nEAS)         
        REAL*8 Weas(6,nEAS) !Weas matric that computes T0(6x6)*Meas(6x30 as maximum)
        REAL*8 EASFor !EAS formulation
        INTEGER EASlocal !paramater
        REAL*8 T0mat(6,6)
!--------------------------------        
        REAL*8 EASstrain(nEAS,1) !enhancing parameters
        REAL*8 Etilde(ntens,1)
        REAL*8 FeEAS(nEAS,1)
        REAL*8 FEAS(nEAS,1)
        REAL*8 KedEAS(ndim*nnode,nEAS)
        REAL*8 KdEAS(ndim*nnode,nEAS)
        REAL*8 KeEASEAS(nEAS,nEAS)
        REAL*8 KEASEAS(nEAS,nEAS)
        REAL*8 KeEASd(nEAS,ndim*nnode)
        REAL*8 KEASd(nEAS,ndim*nnode)
        REAL*8 KeEASEASinv(nEAS,nEAS) !inverse of KeEASEAS
!-------------------------------------------------------
!-----------       ANS SHEAR parameter    --------------
!-------------------------------------------------------
        REAL*8 ANSshearpar      
        REAL*8 gkovr1q1(nANSshear,ndim,ndim)
        REAL*8 gkovc1q1(nANSshear,ndim,ndim)
        REAL*8 gkonr1q1(nANSshear,ndim,ndim)
        REAL*8 gkonc1q1(nANSshear,ndim,ndim)        
        REAL*8 gmkovr1q1(nANSshear,ndim,ndim)
        REAL*8 gmkovc1q1(nANSshear,ndim,ndim)
        REAL*8 gmkonr1q1(nANSshear,ndim,ndim)
        REAL*8 gmkonc1q1(nANSshear,ndim,ndim)        
        REAL*8 detr1(nANSshear)
        REAL*8 detc1(nANSshear)       
        REAL*8 shapef1q1(nANSshear,nnode)
        REAL*8 dshape1q1(nANSshear,ndim,nnode)      
        REAL*8 gkovr2q2(nANSshear,ndim,ndim)
        REAL*8 gkovc2q2(nANSshear,ndim,ndim)
        REAL*8 gkonr2q2(nANSshear,ndim,ndim)
        REAL*8 gkonc2q2(nANSshear,ndim,ndim)        
        REAL*8 gmkovr2q2(nANSshear,ndim,ndim)
        REAL*8 gmkovc2q2(nANSshear,ndim,ndim)
        REAL*8 gmkonr2q2(nANSshear,ndim,ndim)
        REAL*8 gmkonc2q2(nANSshear,ndim,ndim)        
        REAL*8 detr2(nANSshear)
        REAL*8 detc2(nANSshear)        
        REAL*8 shapef2q2(nANSshear,nnode)
        REAL*8 dshape2q2(nANSshear,ndim,nnode)        
        REAL*8 frq(2)
        REAL*8 fsq(2) 
!------------------------------------------------------- 
!-------------------------------------------------------
!-----------       ANS TRAP. parameter    --------------
!-------------------------------------------------------
        REAL*8  ANStrappar       
        REAL*8 gkovrT(nANStrap,ndim,ndim)
        REAL*8 gkovcT(nANStrap,ndim,ndim)
        REAL*8 gkonrT(nANStrap,ndim,ndim)
        REAL*8 gkoncT(nANStrap,ndim,ndim)        
        REAL*8 gmkovrT(nANStrap,ndim,ndim)
        REAL*8 gmkovcT(nANStrap,ndim,ndim)
        REAL*8 gmkonrT(nANStrap,ndim,ndim)
        REAL*8 gmkoncT(nANStrap,ndim,ndim)        
        REAL*8 detrT(nANStrap)
        REAL*8 detcT(nANStrap)        
        REAL*8 shapefT(nANStrap,nnode)
        REAL*8 dshapeT(nANStrap,ndim,nnode)        
        REAL*8 frT(4) 
!------------------------------------------------------- 
!-------------------------------------------------------
      real*8 FeEASn(nEAS,1)
      real*8 KeEASEASn(nEAS,nEAS),KeEASEASinvn(nEAS,nEAS)
      real*8 KeEASdn(nEAS,ndim*nnode) 
      real*8 KedEASn(ndim*nnode,nEAS)
      real*8 Keds(ndim*nnode,nnode),Kesd(nnode,ndim*nnode)
      real*8 KeEASs(nEAS,nnode),KeEASsn(nEAS,nnode)
      real*8 KesEAS(nnode,nEAS),sepf(nnode,nnode),pepf(nnode,1)  
	  REAL*8 Hn
	  REAL*8 coordx !crack length
	  REAL*8 alphT, alphn, alphBn,alphB,alph
	  REAL*8 Fat_deg
      !--------------------------------
      real*8 pf_inc(8,1),pf_v(8)
        real*8 Ls0,Gc,Flag_Fat,Kmodel2,kres
        real*8 EASmodes(nEAS,1)
		real*8 DeltaEASmodes(nEAS,1)
      REAL*8 kGradPhi(ndim,nnode)
      REAL*8 Boppf(ndim,nnode)
      REAL*8 kPhinterp,kPhinterpM,AvgPF          !phase field
      REAL*8 Pf_vn          !phase field at the previous increment
	    REAL*8 DefGrad(3,3),DefGrad_per(3,3),DetF,DefGradInv(3,3) !step n+1
	  REAL*8 DefGradn(3,3),DetFn,DefGradInvn(3,3)  !step n
!--------------------------------


        !matrices
        !Coupling displacement--phase field
        REAL*8 kds(ndim*nnode,nnode)
        REAL*8 ksd(nnode,ndim*nnode)
        !--------------------------------
        !Coupling EAS--phase field
        REAL*8 kEASs(nEAS,nnode)
        REAL*8 ksEAS(nnode,nEAS)
        !--------------------------------
        !Phase field -- phase field
        REAL*8 spf(nnode,nnode)
        REAL*8 ppf(nnode,1)
!-------------------------------------------------------
!-------------------------------------------------------
	     !initialization of residual and tangent in abaqus format
		rhs(:,1) = 0.d0 
		amatrx(:,:) =  0.d0
        !initialization of stiffness matrix and internal force vector
		 Fed(:) = 0.d0
		 Kedd(:,:) = 0.d0
         Kegdd(:,:) = 0.d0
        ! ------- Initialization due to EAS ------------------
            FeEAS(:,1) = 0.d0
            FeEASn(:,1) = 0.d0
        KeEASEAS(:,:) = 0.d0 
        KeEASEASinv(:,:) = 0.d0  
        KeEASEASn(:,:) = 0.d0 
        KeEASEASinvn(:,:) = 0.d0
         KeEASd(:,:) = 0.d0  
        KeEASdn(:,:) = 0.d0 
         KedEAS(:,:) = 0.d0  
        KedEASn(:,:) = 0.d0 
         ! ------- Initialization due to PHASE FIELD ------------------     
			Keds(:,:) = 0.0d0
			Kesd(:,:) = 0.0d0
			KeEASs(:,:) = 0.0d0
			KesEAS(:,:) = 0.0d0
			sepf(:,:) = 0.0d0
			pepf(:,1) = 0.0d0
        !------------------------------------------------------
         !Initialize reference and current coordinates
	   xref(:,:) = 0.d0
       xcur(:,:) = 0.d0  
       !xn(:,:) = 0.d0 !previous converged state	  
        !READING PROCEDURE OF ELEMENT COORDINATES
	    xref(:,:) = coords(:,:)  ! (Xref now has coordniates)
		
! ---------- Begin for update the configuration
! Update the current configuration (Deformation mapping)					 
		do i=1,3
		do j=1,8
	!displacement of node 1-8
	  xcur(i,j) = xref (i,j) + u(i+(j-1)*3)
	 enddo
	  enddo 
    ! Increment of u from i-1 to i (Newton-Raphson)         
			 do i=1,nnode*ndim
			     DeltaU(i,1)=du(i,1)
			 end do      
    ! Increment of s from i-1 to i (Newton-Raphson)
      do i=1,nnode
         pf_v(i)=u(i+nnode*ndim) ! current s 
      enddo    


      !------------------------------------------
      !Reading ELEMENT PROPERTIES FROM INPUT FILE
      !------------------------------------------
      !Formulation due to EAS
      EASpar = props(1) !EAS paramater 1: active 0:innactive
      EASFor = props(2)
      ! ANS 
      ANSshearpar  =props(3) !ANS for transverse shear locking
      ANStrappar =  props(4) !ANS for trapezoidal locking
      !Material paramateres
      emod = props(5)
      enu = props(6)
      
      !PHASE FIELD ENTRIES
	Ls0 = props(7)
	Gc = props(8)
	Flag_Fat = props(9)
	Kmodel2 =  props(10)
	kres= props(11)
	!fatigue
	alphT=Gc/(2.d0*6.d0*Ls0)
      !--------------------------------
      !-------- Compute Guass point information---------------
      call kgaussS10(ngauss,posgp,weightgp)
      call kgaussS10(ngausst,posgpt,weightgpt)
      !--------------------------------
      !--------------------------------------
      !recover Previous EASmodes  
     	do m=1,nEAS
     	  EASmodes(m,1) = 0.d0
		EASmodes(m,1) = svars(m) !recovering the last alpha
	   end do
      
      
          do i=1,nEAS
                FeEASn(i,1)=0.d0
                FeEASn(i,1)=svars(i+7) !This is valid for 7 EAS modes
           end do           
          
		  		  
           do i=1,nEAS !& parameters from 15 to 63 
		   do j=1,7
                KeEASEASinvn(j,i)=0.d0
                KeEASEASinvn(j,i)=svars(i+14+(j-1)*7) ! Just a compact way for notation.
           enddo
		   enddo
		  

           do i=1,ndim*nnode !& parameters from 64 to 231 
		   do j=1,7
                KeEASdn(j,i)=0.d0
                KeEASdn(j,i)=svars(i+63+(j-1)*24)
           enddo
		   enddo
 !For fatigue       
           Hn=0.d0
           Hn=svars(232)
           alphn=svars(233)
		   alphBn=svars(234)
    
           
           !In total we have 234 svars
           !----  
           ! EASpar (7)            1-7
           ! FeEAS  (7)            8-14
           ! KeEASEASinv (7x7=49) 15-63
           ! KeEASd  (7x24=168)   64-231
		   ! Hn                  232
		   ! Alphn                  233
		   ! AlphBn                  234
           
      !2 Compute DeltaEASmodes   
      call kDeltaEASvect(DeltaEASmodes,FeEASn,KeEASdn,KeEASEASinvn,
     * deltau,ndim,nnode,nEAS)
     
       do i=1,nEAS
                !EASmodes(i,1)=0.0d0
                EASmodes(i,1)=EASmodes(i,1)+DeltaEASmodes(i,1)
       end do
        !END of update the internal values due to EAS at the previous converged  configuration
      !---------------------------
      !-- Start EAS Now -------------------
      !---------------------------
	if (EASpar.ne.0.d0) then

		!Evaluation at the element center
      call kShapeFunctS10(0.d0,0.d0,0.d0,ndim,nnode,shapef0,dshape0)

      !METRICS AT THE ELEMENT center for EAS: reference and current configuration
       call ks8metShellBody10(xref,shapef0,dshape0,0.d0,nnode,ndim,
     *   gkovr0,gkonr0,gmkovr0,gmkonr0,detr0)

!      call ks8metShellBody10(xcur,shapef0,dshape0,0.d0,nnode,ndim,
!     *   gkovc0,gkonc0,gmkovc0,gmkonc0,detc0)  

	end if
	!---------------------------
      !--- END EAS ---------------
      !---------------------------
      !---------------------------
      !--  ANS transv. shear -----
      !---------------------------
      if (ANSshearpar.ne.0.d0) then

      call kANScolpointsShear(gkovr1q1,gkovc1q1,gkonr1q1,gkonc1q1,
     * gmkovr1q1,gmkovc1q1,gmkonr1q1,gmkonc1q1,detr1,detc1,shapef1q1,
     * dshape1q1,gkovr2q2,gkovc2q2,gkonr2q2,gkonc2q2,gmkovr2q2,
     * gmkovc2q2,gmkonr2q2,gmkonc2q2,detr2,detc2,shapef2q2,dshape2q2,
     * xref,xcur,ndim,nnode,nANSshear)

      endif
      !---------------------------
      !-- END ANS transv. shear --
      !---------------------------

      !---------------------------
      !--  ANS trapez.   ---------
      !---------------------------
      if (ANStrappar.ne.0.d0) then

      call  kANScolpointsTrap(gkovrT,gkovcT,gkonrT,gkoncT,gmkovrT,
     * gmkovcT,gmkonrT,gmkoncT,detrT,detcT,shapefT,dshapeT,
     * xref,xcur,ndim,nnode,nANStrap)

      endif
      !---------------------------
      !-- END ANS trapez. --------
      !---------------------------


      kinkt = 0 !counter
      !LOOP OVER THE IP
       !Loop over the eta-direction
	  kPhinterpM=0.d0
      do i=1,ngauss
        eta =  posgp(i)
        weta =  weightgp(i)

        !Loop over the phi-direction
        do j=1,ngauss
            phi =  posgp(j)
            wphi =  weightgp(j)

            !Loop over the thickness
            do k=1,ngausst
                xi = posgpt(k)
                wxi =  weightgpt(k)
                kinkt = kinkt+ 1

!                !Maximum phase field
                Pf_vn=0.d0
                !Pf_vn =  svars(234 + kinkt)
				
			isvinc= (kinkt-1)*nsvarsip !nsvarsip is the number of state dependent variables
                
                do l = 1,nsvarsip !for each integration points if the element
				  !print*,l+isvinc,"l+isvinc"
                  statevLocal(l)=svars(234+l+isvinc) !SVDs for the nodes
                end do
                
                !Psin=statevLocal(1)
                !Hn=statevLocal(2)
                !alphn=statevLocal(3)
                !alphBn=statevLocal(4)
                Pf_vn=statevLocal(5)

       !Compute shape functions and derivatives
           call kShapeFunctS10(eta,phi,xi,ndim,nnode,shapef,dshape)


      !Compute metrics at the reference and current configurations
           call ks8metShellBody10(xref,shapef,dshape,xi,nnode,ndim,
     *     gkovr,gkonr,gmkovr,gmkonr,detr) !reference


           call ks8metShellBody10(xcur,shapef,dshape,xi,nnode,ndim,
     *     gkovc,gkonc,gmkovc,gmkonc,detc)!current
 
      !---------------------------------------
      !Compute B-operator: DISPLACEMENT DERIVATIVE
           call kBoperatorS10(Bop,gkovc,dshape,nnode,ndim)
		         !------------------------------------

      Jacobtr=weta*wphi*wxi*detr; !compute the constant to multiply for integration
	  ! compute the deformation gradient in cartesian
	  call kDefGradient(DefGrad,DetF,DefGradInv,gkovc,gkonr)
		
		
		
      !Constitutive Law
      call kMatLawS10(Cmat,gmkovc,gmkovr,gmkonr,gkonr, gkovr,
     *   DefGrad,DetF,DefGradInv, emod,enu,stress,strain, 
     * ntens,KFreeEner,kFreeEner_pos,kFreeEner_neg)
      !Material stiffness matrix
      call kStiffMatrixS10(Kdd,Fd,nnode,ndim,ndof,Jacobtr,
     *     Cmat,stress,Bop,ntens)
      !---------------------------------------
      !------------------------------------
      !--------  ANS shear comp. ----------
      !------------------------------------
	  if (Pf_vn.lt.PFTol) then
	  if (ANSshearpar.ne.0.d0) then

        call kshapefANSshear(eta,phi,frq,fsq)

        call kANSshearS10(Bop,ndim,nnode,ntens,ndof,gkovc1q1,
     *  gkovc2q2,frq,fsq,dshape1q1,dshape2q2,nANSshear)
      endif

      !------------------------------------
      !------  ANS trapez. comp. ----------
      !------------------------------------
      if (ANStrappar.ne.0.d0) then
        call kshapefANStrap(eta,phi,frT)

        call  kANStrapS10(Bop,ndim,nnode,ntens,ndof,gkovcT,
     *  frT,dshapeT,nANStrap)
      endif
	  endif
      !------------------------------------
      !--------  EAS computation ----------
      !------------------------------------
        if (EASpar.ne.0.d0) then

       !Interpolation of the EAS paramaters
       call kEASMlocal(eta,phi,xi,EASfor,Weasloc,EASlocal,nEAS)
       !Transformation from the element center coordinates
       call kEASTransform(Weas,T0mat,gkovr0,gkonr0,Weasloc,EASlocal,
     * ndim,nnode,detr0,detr,nEAS)
       call kEtilde(Etilde,Weas,EASmodes,nEAS)
         end if



      !Geometrical stiffness matrix in Ramm way
	  if (Pf_vn.lt.PFTol) then
      if ((ANSshearpar.ne.0.d0).OR.(ANStrappar.ne.0.d0)) then

      call kStiffGeomS10AltANS(Kgdd,shapef,dshape,stress,ndim,
     * nnode,ntens,gkovc1q1,gkovc2q2,frq,fsq,dshape1q1,dshape2q2,
     * nANSshear,ANSshearpar,gkovcT,frT,dshapeT,nANStrap,ANStrappar)

      endif
	  else
      call kStiffGeomS10Alt(Kgdd,shapef,dshape,stress,ndim,nnode,
     *  ntens)	  
	  endif

      !---------------------------------------------------------
      !----------     PHASE FIELD block     --------------------
      !---------------------------------------------------------

      !Compute material gradient of shape functions
      call kGradShapeFunctionS10(kGradPhi,Boppf,
     * dshape,gkonr,nnode,ndim)

      !Compute Phase Field and Gradient --> s and Grad(s)
      call kPhaseInterpS10(kPhinterp,kGradPhi,shapef,Boppf,
     * ndim,nnode,Pf_v,kFreeEner_pos,kFreeEner_neg)
	 
!	    svars(234+kinkt)=kPhinterp
		

	    kPhinterpM=kPhinterpM+kPhinterp

      !write(7,*) "After kPhaseInterpS10"
!     update crack length
       coordx=0.d0
       if (kPhinterp.ge.0.975d0) then
        do ii=1,nnode
         coordx=coordx+shapef(ii)*coords(2,ii)
        enddo
        if (coordx.gt.Ac) then
         Ac=coordx
        endif
       endif

       !Maximum phase field: update of the history variables
       alph=kFreeEner_pos*((1-kPhinterp)**2)   !maybe here
!     Update fatigue history variable
       if ((alph.ge.alphn).and.(dtime.gt.0.d0)) then
        alphB = alphBn+abs(alph-alphn)
       else
        alphB=alphBn
       endif

      !Compute  KKT conditions
	  if (Hn .le. kFreeEner_pos) then
	      Hn = kFreeEner_pos
	  endif
	  !fatigue Degradration
	  if (Flag_Fat.eq.1) then 
	  
	  !fatigue Degradration
       if (alphB.lt.alphT) then
        Fat_deg= 1.d0
       else
        Fat_deg=(2.d0*alphT/(alphB+alphT))**kmodel2
      	 endif
	  
	  endif 
	  
	 
	  	  if (Flag_Fat.eq.2) then 
		!write(7,*)'Fatigue Degradation-2'
	  !fatigue Degradration
       if (alphB.lt.alphT) then
        Fat_deg= 1.d0
       endif 
	   if ((alphT.le.alphB).and.(alphB.lt.alphT*(10)**(1/kmodel2))) then
        Fat_deg=(1-kmodel2* log(alphB/alphT))**2.d0
      endif
	   if(alphB.gt.alphT*(10)**(1/kmodel2)) then 
		  Fat_deg=0.0d0
	  endif 
	  
	  endif
	  
	  
	  
        call kPhaseFieldResS10(ppf,kPhinterp,kGradPhi,Hn,Gc,Ls0,
     *  shapef,Boppf,nnode,ndim,Jacobtr,Fat_deg)

        do ik=1,nnode
             pepf(ik,1) = pepf(ik,1) + ppf(ik,1)
        enddo

       !Compute phase field matrix spf
        call kPhaseFieldStiffS10(spf,Hn,Gc,Ls0,
     *  shapef,Boppf,nnode,ndim,Jacobtr,Fat_deg)

        do ik=1,nnode
           do jk=1,nnode
             sepf(ik,jk) = sepf(ik,jk) + spf(ik,jk)
           enddo
        enddo
!------------------------------------------------------------------------
       do ik=1,nnode*ndim
          Fed(ik) = Fed(ik)  + Fd(ik)*((1.d0-Pf_vn)**2.0d0)
          do jk=1,nnode*ndim
             Kedd(ik,jk) = Kedd(ik,jk)+Kdd(ik,jk)*((1.d0-Pf_vn)**2.0d0)
             Kegdd(ik,jk) = Kegdd(ik,jk)+Kgdd(ik,jk)*Jacobtr*((1.d0-Pf_vn)**2.0d0)
          enddo
       enddo

      !------------------------------------
      !--------  EAS computation:  --------
      !        Residual and Matrices ------
      !------------------------------------
       if (EASpar.ne.0.d0) then

       call kresidualEAS(FEAS,Weas,stress,nEAS)
       do ik=1,nEAS
            FeEAS(ik,1)=FeEAS(ik,1)+FEAS(ik,1)*Jacobtr
       enddo
	   
       !Coupling EAS-displacements
       call kmatricesEAS(KdEAS,KEASEAS,KEASd,Weas,Cmat,Bop,nEAS)
       
       do ik=1,nEAS
          do jk=1,nEAS
         KeEASEAS(ik,jk) =KeEASEAS(ik,jk)+KEASEAS(ik,jk)*Jacobtr
          enddo
       enddo	   
 
       do ik=1,nEAS
          do jk=1,ndim*nnode
         KeEASd(ik,jk) =KeEASd(ik,jk)+KEASd(ik,jk)*Jacobtr
          enddo
       enddo	   

       do ik=1,ndim*nnode
          do jk=1,nEAS
        KedEAS(ik,jk) = KedEAS(ik,jk) + KdEAS(ik,jk)*Jacobtr
          enddo
       enddo  

       endif !EAS
       !----------------------------
	       statevLocal(1)= kFreeEner_pos
		  statevLocal(2)= strain(4)
		  statevLocal(3)= Ac
		  statevLocal(4)= stress(4)
		  statevLocal(5)= kPhinterp
            
            do l = 1,nsvarsip !for each integration points if the element
              !print*,l+isvinc,"l+isvinc"
              svars(234+isvinc+l)=statevLocal(l) !SVDs for the nodes
              USRVAR(jelem,l,kinkt)=statevLocal(l)
            end do
	 

            enddo ! end loop over xi
        enddo ! end loop over phi
      enddo ! end loop over eta


      !--- Geometrically nonlinear
      NLGEOMflag  = 1.d0
      if (NLGEOMflag.ne.0.d0) then

      !Summ material and geometric contributions
       do i=1,ndim*nnode
            do j=1,ndim*nnode
               Kedd(i,j) = Kedd(i,j)  + Kegdd(i,j)
            enddo
       enddo
       
      endif
       !------


      !------------------------------------
      !--------  EAS computation:  --------
      !---  Static condensation process ---
      !------------------------------------
      AvgPF=0.d0
	  AvgPF=kPhinterpM/nnode
!	  if (SintA .le. 0.2d0) then
	  if (EASpar.ne.0.d0) then
	  
      if (AvgPF .lt. PFTol) then
	  
      call InverMatrixS10(KeEASEAS,KeEASEASinv,nEAS)
        call kStaticCondenS10(Kedd,Fed,KedEAS,KeEASEASinv,
     * KeEASd,FeEAS,ndim,nnode,nEAS)

       endif
	   endif

       !-------------------------------
       !--------------------------------
       !-------------------------------
			rhs(1:24,1)  = rhs(1:24,1) - Fed(1:24)
		 amatrx(1:24,1:24) = amatrx(1:24,1:24) + Kedd(1:24,1:24)
		 
		 
		 rhs(25:32,1)  = rhs(25:32,1) - pepf(1:8,1)
		 amatrx(25:32,25:32) = amatrx(25:32,25:32) + sepf(1:8,1:8)
		   

  
  !Previous EASmodes 
     	do m=1,nEAS
		svars(m) = EASmodes(m,1) !recovering the last alpha
	    end do
      
           !----------
           do i=1,nEAS
                svars(i+7)=FeEAS(i,1) !This is valid for 3 EAS modes
           end do   
           !----------
		   
		   
           do i=1,nEAS
		   do j=1,7
                svars(i+14+(j-1)*7) =KeEASEASinv(j,i)
           enddo
		   enddo

           do i=1,ndim*nnode  
			do j=1, 7
                svars(i+63+(j-1)*24)=KeEASd(j,i)
           enddo  
           enddo		  

!---------------------------------------------------------
       !historical variable
          svars(232)=Hn
       !fatigue
          svars(233)=alph
          svars(234)=alphB		  
   
      return 
      end

       
      subroutine kgaussS10(ngauss,posgp,weightgp)
      
         !IMPLICIT NONE 
      
      INTEGER ngauss  !number of IP in-plane
      INTEGER ngausst !number of IP over thickness
      
      REAL*8 posgp(ngauss),weightgp(ngauss)
      
      if (ngauss.eq.2) then
      
      	posgp(1) = -0.577350269189626d0
		weightgp(1) = 1.d0
		
		posgp(2) = 0.577350269189626d0
		weightgp(2) = 1.d0
		
      endif
      
      return
      end  
!---------------------------------------------------------
         subroutine  kDefGradient(F,JacobF,Finv,gkovc,gkonr)
           implicit none 
	   
       !Subroutine to compute the shape functions and ther derivatives  
        INTEGER i,j,k,l
        REAL*8 F(3,3),JacobF,Finv(3,3),gkovc(3,3),gkonr(3,3)
        
        !initialization
		F(:,:)=0.0d0
		Finv(:,:)=0.0d0
		JacobF = 0.d0
		

      
        do i=1,3
           do j=1,3
            do k=1,3
            F(j,i) =  F(j,i) + gkovc(i,k)*gkonr(j,k);
            enddo
          enddo
        enddo
	  
	  
	  JacobF = F(1,1)*F(2,2)*F(3,3)+ 
     *  F(1,3)*F(2,1)*F(3,2)+ F(3,1)*F(1,2)*F(2,3) -  
     *      F(3,1)*F(2,2)*F(1,3)-F(3,3)*F(1,2)*F(2,1)-  
     *      F(1,1)*F(2,3)*F(3,2)
     
     
         if(JacobF.gt.0.0d0)  then
	      Finv(1,1)=  (F(2,2) *F(3,3) - F(3,2) *F(2,3))/JacobF
	      Finv(1,2)= -(F(1,2) *F(3,3) - F(1,3) *F(3,2))/JacobF
          Finv(1,3)=  (F(1,2) *F(2,3) - F(2,2) *F(1,3))/JacobF
	      Finv(2,1)= -(F(2,1) *F(3,3) - F(3,1) *F(2,3))/JacobF
	      Finv(2,2)=  (F(1,1) *F(3,3) - F(1,3) *F(3,1))/JacobF
	      Finv(2,3)= -(F(1,1) *F(2,3) - F(2,1) *F(1,3))/JacobF
	      Finv(3,1)=  (F(2,1) *F(3,2) - F(3,1) *F(2,2))/JacobF
	      Finv(3,2)= -(F(1,1) *F(3,2) - F(3,1) *F(1,2))/JacobF
	      Finv(3,3)=  (F(1,1) *F(2,2) - F(1,2) *F(2,1))/JacobF
	      endif   
		 return
         end
!--------------------------------------------------------------		 
!---------------------------------------------------------   
       subroutine kShapeFunctS10(r,s,t,ndim,nnode,shapef,dshape)
       
	      !IMPLICIT NONE 
	   
       !Subroutine to compute the shape functions and ther derivatives  
        INTEGER nnode, ndim
        REAL*8 shapef(nnode),r,s,t, dshape(ndim,nnode)
        
        !initialization
		shapef(:)=0.0d0
		dshape(:,:)=0.0d0

      shapef(1) = ((1-R)*(1+S)*(1-T))/8.0d0
      dshape(1,1) =  (-(1+S)*(1-T))/8.0d0
      dshape(2,1) =  ((1-R)*(1-T))/8.0d0
      dshape(3,1) =  (-(1-R)*(1+S))/8.0d0
 
	  shapef(2) = ((1-R)*(1-S)*(1-T))/8.0d0
	  dshape(1,2) = (-(1-S)*(1-T))/8.0d0
	  dshape(2,2) = (-(1-R)*(1-T))/8.0d0
	  dshape(3,2) = (-(1-R)*(1-S))/8.0d0
    
	  shapef(3) = ((1+R)*(1-S)*(1-T))/8.0d0
      dshape(1,3) = ((1-S)*(1-T))/8.0d0
      dshape(2,3) = (-(1+R)*(1-T))/8.0d0
      dshape(3,3) = (-(1+R)*(1-S))/8.0d0
    
	  shapef(4) = ((1+R)*(1+S)*(1-T))/8.0d0
      dshape(1,4) =  ((1+S)*(1-T))/8.0d0
      dshape(2,4) =  ((1+R)*(1-T))/8.0d0
      dshape(3,4) = (-(1+R)*(1+S))/8.0d0

      shapef(5) = ((1-R)*(1+S)*(1+T))/8.0d0
      dshape(1,5) = (-(1+S)*(1+T))/8.0d0
      dshape(2,5) = ((1-R)*(1+T))/8.0d0
      dshape(3,5) = ((1-R)*(1+S))/8.0d0

	  shapef(6) = ((1-R)*(1-S)*(1+T))/8.0d0
	  dshape(1,6) = (-(1-S)*(1+T))/8.0d0
	  dshape(2,6) = (-(1-R)*(1+T))/8.0d0
	  dshape(3,6) = ((1-R)*(1-S))/8.0d0

	  shapef(7) = ((1+R)*(1-S)*(1+T))/8.0d0
      dshape(1,7) = ((1-S)*(1+T))/8.0d0
      dshape(2,7) = (-(1+R)*(1+T))/8.0d0
      dshape(3,7) = ((1+R)*(1-S))/8.0d0
	  
	  shapef(8) = ((1+R)*(1+S)*(1+T))/8.0d0
      dshape(1,8) = ((1+S)*(1+T))/8.0d0
      dshape(2,8) = ((1+R)*(1+T))/8.0d0
      dshape(3,8) = ((1+R)*(1+S))/8.0d0
     
       return
       end 
            
!---------------------------------------------------------    
!---------------------------------------------------------         
!--------------------------------------------------------- 
!--------------------------------------------------------- 
      subroutine kDeltaEASvect(DeltaEAS,FeEASn,KeEASdn,KeEASEASinvn,
     * DeltaU,ndim,nnode,nEAS)
     
         !IMPLICIT NONE
      
      REAL*8 FeEASauxn(nEAS,1), DeltaEAS(nEAS,1), FeEASn(nEAS,1),DeltaU(ndim*nnode,1)   
      REAL*8 KeEASdn(nEAS,ndim*nnode), KeEASEASinvn(nEAS,nEAS)
      INTEGER ndim,nnode,nEAS, i 
      
      !initialization
      FeEASauxn(:,1)=0.d0
      DeltaEAS(:,1)=0.d0
      call kmatrixmultiS10(FeEASauxn,nEAS,1,KeEASdn,nEAS,ndim*nnode,
     * DeltaU,ndim*nnode,1)
      
      do i=1,nEAS
      FeEASauxn(i,1) = FeEASauxn(i,1) + FeEASn(i,1)
      enddo

      call kmatrixmultiS10(DeltaEAS,nEAS,1,KeEASEASinvn,nEAS,nEAS,
     * FeEASauxn,nEAS,1)

      do i=1,nEAS
      DeltaEAS(i,1) = -1.d0*DeltaEAS(i,1)
      enddo

      RETURN
      END
!---------------------------------------------------------    
!---------------------------------------------------------
       subroutine ks8metShellBody10(x,shapef,dshape,xi,nnode,ndim,
     *   gkov,gkon,gmkov,gmkon,detgkon) 
  
        INTEGER nnode, ndim
        REAL*8 xi, x(ndim,nnode), betr, shapef(nnode),dshape(ndim,nnode)
        REAL*8 gkov(3,3),gkon(3,3), gmkov(3,3),gmkon(3,3),detr,detgkon,scalar
        REAL*8 detgmkov, t(3,3), transcart(6,6), tof(3,3)
        REAL*8 gkoni(3,3), contra(3,3), covariant(3,3), cro(3)
  
        
        !initialization
              gkov(:,:) = 0.0d0  
              gkon(:,:) = 0.0d0  
              gmkov(:,:) = 0.0d0  
              gmkon(:,:) = 0.0d0 
              gkoni(:,:) = 0.0d0 
              gkoni(:,:) = 0.0d0 
			contra(:,:)=0.0d0
			covariant(:,:)=0.0d0		
          !write(7,*) "x in metrics", x

      do i=1,ndim
         do j =1,ndim
           do k=1,nnode     
       gkov(j,i) = gkov(j,i) + dshape(i,k)*x(j,k)
             enddo
         enddo
      enddo  
			!gkov= Covariant basis vector
      covariant(:,:)=gkov(:,:)
      do i=1,ndim
        do j=1,ndim
      gkon(i,j) = gkov(i,j)  
        enddo
      enddo
      

      detgkon=gkon(1,1)*gkon(2,2)*gkon(3,3)+
     *        gkon(1,3)*gkon(2,1)*gkon(3,2)+
     *        gkon(3,1)*gkon(1,2)*gkon(2,3)- 
     *        gkon(3,1)*gkon(2,2)*gkon(1,3)-
     *        gkon(3,3)*gkon(1,2)*gkon(2,1)-
     *        gkon(1,1)*gkon(2,3)*gkon(3,2)
     
     
      	if(detgkon.ne.0.d0) then
      gkoni(1,1)=  (gkon(2,2) *gkon(3,3) - gkon(3,2) *gkon(2,3))/detgkon
      gkoni(1,2)= -(gkon(1,2) *gkon(3,3) - gkon(1,3) *gkon(3,2))/detgkon
      gkoni(1,3)=  (gkon(1,2) *gkon(2,3) - gkon(2,2) *gkon(1,3))/detgkon
      gkoni(2,1)= -(gkon(2,1) *gkon(3,3) - gkon(3,1) *gkon(2,3))/detgkon
      gkoni(2,2)=  (gkon(1,1) *gkon(3,3) - gkon(1,3) *gkon(3,1))/detgkon
      gkoni(2,3)= -(gkon(1,1) *gkon(2,3) - gkon(2,1) *gkon(1,3))/detgkon
      gkoni(3,1)=  (gkon(2,1) *gkon(3,2) - gkon(3,1) *gkon(2,2))/detgkon
      gkoni(3,2)= -(gkon(1,1) *gkon(3,2) - gkon(3,1) *gkon(1,2))/detgkon
      gkoni(3,3)=  (gkon(1,1) *gkon(2,2) - gkon(1,2) *gkon(2,1))/detgkon
      else 
       !!write(7,*)'Warning: element  has negative Jacobtr in detgkon'   
      end if
      
      do i=1,ndim
         do j=1,ndim 
        gkon (i,j) = gkoni(i,j)   
          end do 
        end do

c  calculating the transponse

       do i=1,ndim
          do j=i+1,ndim         
           scalar = gkon(j,i)
           gkon (j,i) =gkon(i,j)
           gkon(i,j)=scalar  
          end do 
      end do

c ------- covariant metric tensor -------
        do i=1,ndim
                do j=1,ndim
                        do k=1,ndim 
            gmkov(i,j) =gmkov(i,j) + gkov(k,i)*gkov(k,j)
                        end do
                end do
        end do
c  -------   making it symmetric  -------
        gmkov(2,1) = gmkov(1,2)
        gmkov(3,1) = gmkov(1,3)
        gmkov(3,2) = gmkov(2,3)
c ------- calculating the inverse gmkon  -------
      detgmkov=gmkov(1,1)*gmkov(2,2)*gmkov(3,3)+
     * gmkov(1,3)*gmkov(2,1)*gmkov(3,2)+
     * gmkov(3,1)*gmkov(1,2)*gmkov(2,3)-
     * gmkov(3,1)*gmkov(2,2)*gmkov(1,3)-
     * gmkov(3,3)*gmkov(1,2)*gmkov(2,1)-
     * gmkov(1,1)*gmkov(2,3)*gmkov(3,2)
	 
c *** Now the inverse **
         if(detgmkov.gt.0.d0) then
      
      gmkon(1,1) = (gmkov(2,2) *gmkov(3,3) -
     *              gmkov(3,2) *gmkov(2,3))/detgmkov
     
      gmkon(1,2) = -(gmkov(1,2) *gmkov(3,3) - 
     *               gmkov(1,3) *gmkov(3,2))/detgmkov
     
      gmkon(1,3) =  (gmkov(1,2) *gmkov(2,3) - 
     *               gmkov(2,2) *gmkov(1,3))/detgmkov
     
      gmkon(2,1) = -(gmkov(2,1) *gmkov(3,3) -
     *               gmkov(3,1) *gmkov(2,3))/detgmkov
     
      gmkon(2,2) =  (gmkov(1,1) *gmkov(3,3) -
     *               gmkov(1,3) *gmkov(3,1))/detgmkov
     
      gmkon(2,3) = -(gmkov(1,1)*gmkov(2,3) - 
     *               gmkov(2,1)*gmkov(1,3))/detgmkov
     
      gmkon(3,1) =  (gmkov(2,1)*gmkov(3,2) - 
     *               gmkov(3,1)*gmkov(2,2))/detgmkov
     
      gmkon(3,2) = -(gmkov(1,1) *gmkov(3,2) - 
     *               gmkov(3,1) *gmkov(1,2))/detgmkov
     
      gmkon(3,3) =  (gmkov(1,1) *gmkov(2,2) - 
     *               gmkov(1,2) *gmkov(2,1))/detgmkov
        else 
         !!write(7,*)'Warning: element has negative Jacobtr in detgmkov' 
        end if     


  
        return
        end
!---------------------------------------------------------    
!---------------------------------------------------------
!---------------------------------------------------------
      subroutine kBoperatorS10(Bop,gkov,dshape,nnode,ndim)

      INTEGER ndim,nnode,node_st,inode
      REAL*8 gkov(3,3), Bop_curv(6, ndim*nnode),dshape(ndim,nnode)
      REAL*8 Bop(6,ndim*nnode), transcart(6,6),Bop_tra(6,ndim*nnode)
	  
	  Bop(:,:)=0.d0        
      node_st=1
	do inode=1,nnode
	
      Bop(1,node_st)=       dshape(1,inode)*gkov(1,1)
      Bop(1,node_st+1)=     dshape(1,inode)*gkov(2,1)
      Bop(1,node_st+2)=     dshape(1,inode)*gkov(3,1)
       !----------------------------------------
      Bop(4,node_st)=      dshape(2,inode)*gkov(1,2)
      Bop(4,node_st+1)=    dshape(2,inode)*gkov(2,2)	
      Bop(4,node_st+2)=    dshape(2,inode)*gkov(3,2)
       !----------------------------------------
      Bop(6,node_st) =    dshape(3,inode)*gkov(1,3) 
      Bop(6,node_st+1) =  dshape(3,inode)*gkov(2,3) 
      Bop(6,node_st+2) =  dshape(3,inode)*gkov(3,3) 
      !----------------------------------------
      Bop(2,node_st) =    dshape(1,inode)*gkov(1,2) + 
     *                    dshape(2,inode)*gkov(1,1)
     
      Bop(2,node_st+1) =  dshape(1,inode)*gkov(2,2) +
     *                    dshape(2,inode)*gkov(2,1)
     
      Bop(2,node_st+2) =  dshape(1,inode)*gkov(3,2) +  
     *                    dshape(2,inode)*gkov(3,1) 
      !----------------------------------------
      Bop(3,node_st) =    dshape(1,inode)*gkov(1,3) + 
     *                    dshape(3,inode)*gkov(1,1)
     
      Bop(3,node_st+1) =  dshape(1,inode)*gkov(2,3) +
     *                    dshape(3,inode)*gkov(2,1)
     
      Bop(3,node_st+2) =  dshape(1,inode)*gkov(3,3) + 
     *                    dshape(3,inode)*gkov(3,1)
      !----------------------------------------
      Bop(5,node_st )  =   dshape(2,inode)*gkov(1,3) + 
     *                     dshape(3,inode)*gkov(1,2)
     
      Bop(5,node_st+1) =   dshape(2,inode)*gkov(2,3) + 
     *                      dshape(3,inode)*gkov(2,2)
     
      Bop(5,node_st+2) =   dshape(2,inode)*gkov(3,3) + 
     *                     dshape(3,inode)*gkov(3,2)
	node_st=node_st+ndim		
      enddo

      return
      end
!---------------------------------------------------------    
!--------------------------------------------------------- 
!---------------------------------------------------------       
!---------------------------------------------------------        
      subroutine kANScolpointsTrap(gkovrT,gkovcT,gkonrT,gkoncT,gmkovrT,
     * gmkovcT,gmkonrT,gmkoncT,detrT,detcT,shapefT,dshapeT,
     * xref,xcur,ndim,nnode,nANStrap)  
       
        INTEGER i,j,iq, ndim,nnode,nANStrap
                
        REAL*8 gkovrT(nANStrap,ndim,ndim), gkovcT(nANStrap,ndim,ndim)
        REAL*8 gkonrT(nANStrap,ndim,ndim), gkoncT(nANStrap,ndim,ndim)
        REAL*8 gmkovcT(nANStrap,ndim,ndim), gmkovrT(nANStrap,ndim,ndim)
        REAL*8 gmkonrT(nANStrap,ndim,ndim),gmkoncT(nANStrap,ndim,ndim)        
        REAL*8 detrT(nANStrap), detcT(nANStrap), shapefT(nANStrap,nnode)
        REAL*8 dshapeT(nANStrap,ndim,nnode), xref(ndim,nnode),xcur(ndim,nnode)  
        REAL*8 xrT(nANStrap),xsT(nANStrap), shapefTi(nnode), dshapeTi(ndim,nnode)
        REAL*8 gkovrTi(ndim,ndim),gkonrTi(ndim,ndim),gmkovrTi(ndim,ndim) 
        REAL*8 gmkonrTi(ndim,ndim), gkovcTi(ndim,ndim),gkoncTi(ndim,ndim) 
        REAL*8 gmkovcTi(ndim,ndim),gmkoncTi(ndim,ndim), detrTi,detcTi
		REAL*8 contrai(3,3), transcarti(6,6),contraci(3,3), transcartci(6,6)
        
        !initialization 

        xrT(:) = 0.d0
       xsT(:) = 0.d0
       detrT(:) = 0.d0
       detcT(:) = 0.d0     
        gkovrT(:,:,:) = 0.d0
        gkovcT(:,:,:) = 0.d0
        gkonrT(:,:,:) = 0.d0
        gkoncT(:,:,:) = 0.d0
        gmkovrT(:,:,:) = 0.d0
        gmkovcT(:,:,:) = 0.d0
        gmkonrT(:,:,:) = 0.d0
        gmkoncT(:,:,:) = 0.d0   
        shapefT(:,:) = 0.d0    
        dshapeT(:,:,:) = 0.d0        

        !-----------------------
        xrT(1) =  1.d0
        xrT(2) = -1.d0  

        xsT(1) = 1.d0
        xsT(2) = 1.d0 
        
        xrT(3) = -1.d0
        xrT(4) =  1.d0  

        xsT(3) = -1.d0
        xsT(4) = -1.d0 
        !-----------------------
      
      do  iq =1,nANStrap
      !Compute shape functions and derivatives
      call kShapeFunctS10(xrT(iq),xsT(iq),0.d0,ndim,nnode,
     * shapefTi,dshapeTi) 
    
        do i=1,nnode
        shapefT(iq,i) = shapefTi(i)
        enddo 
        
       do i=1,ndim
        do j=1,nnode
        dshapeT(iq,i,j) =  dshapeTi(i,j) 
         enddo
       enddo
      
      !Compute metrics reference 
      call ks8metShellBody10(xref,shapefTi,dshapeTi,0.d0,nnode,ndim,
     * gkovrTi,gkonrTi,gmkovrTi,gmkonrTi,detrTi)
    
      !Compute metrics current
      call ks8metShellBody10(xcur,shapefTi,dshapeTi,0.d0,nnode,ndim,
     * gkovcTi,gkoncTi,gmkovcTi,gmkoncTi,detcTi)
   
      do i=1,ndim
        do j=1,ndim
        gkovrT(iq,i,j) =   gkovrTi(i,j)   
        gkovcT(iq,i,j) =   gkovcTi(i,j) 
    
        gkonrT(iq,i,j) =   gkonrTi(i,j)
        gkoncT(iq,i,j) =   gkoncTi(i,j)
    
        gmkovrT(iq,i,j) = gmkovrTi(i,j)
        gmkovcT(iq,i,j) = gmkovcTi(i,j)
    
        gmkonrT(iq,i,j) = gmkonrTi(i,j)
        gmkoncT(iq,i,j) = gmkoncTi(i,j)
    
        detrT(iq) = detrTi
        detcT(iq) = detcTi
        enddo
      enddo
      
      
      enddo
      
      RETURN
      END 
!---------------------------------------------------------
!---------------------------------------------------------     
!---------------------------------------------------------    
!--------------------------------------------------------- 
      subroutine kEASMlocal(r,s,t,EASfor,Weasloc,EASlocal,nEAS)
         
         !IMPLICIT NONE      
       
      REAL*8  r,s,t, EASFor, Weasloc(6,nEAS)    
      INTEGER nEAS,indexcol, EASlocal 
      INTEGER i,j, rr,rs,rt,ss,st,tt
       
       !initialization 
      Weasloc(:,:) = 0.d0                    
       rr = 1
       rs = 2
       rt = 3
       ss = 4
       st = 5 
       tt = 6
       
       indexcol = 1
      
      Weasloc(rr,1) = r 
      Weasloc(rr,2) = r*s
      Weasloc(ss,3) = s
      Weasloc(ss,4) = r*s
      Weasloc(tt,5) = t 
      Weasloc(tt,6) = t*r 
      Weasloc(tt,7) = t*s

	indexcol = indexcol + 7 
         
         RETURN
         END  
!---------------------------------------------------------    
!--------------------------------------------------------- 
        
      subroutine kEASTransform(Weas,T0mat,gkovr,gkonr0,Weasloc,EASlocal,
     *  ndim,nnode,detr0,detr,nEAS)
      
         !IMPLICIT NONE

      REAL*8 gkovr(ndim,ndim), gkonr0(ndim,ndim),T0mat(6,6),transcart(6,6)
      REAL*8 Weasloc(6,nEAS),Weas(6,nEAS), Weas_cart(6,nEAS)
      REAL*8 factor0,detr0,detr, TOL     
      REAL*8 t11,t12,t13,t21,t22,t23,t31,t32,t33
      INTEGER ndim,nnode,EASlocal,nEAS,i,j,k,l
 
      TOL = 1.0d-15
        T0mat(:,:)=0.d0
        Weas(:,:) = 0.d0

      factor0= detr0/detr !relation between both determinants

      t11 = 0.d0
      t12 = 0.d0
      t13 = 0.d0
	  
      t21 = 0.d0
      t22 = 0.d0
      t23 = 0.d0

      t31 = 0.d0
      t32 = 0.d0
      t33 = 0.d0
      
      
      do i=1,3
      t11 = t11 + gkonr0(i,1)* gkovr(i,1)
      t12 = t12 + gkonr0(i,1)* gkovr(i,2)
      t13 = t13 + gkonr0(i,1)* gkovr(i,3)

      t21 = t21 + gkonr0(i,2)* gkovr(i,1)
      t22 = t22 + gkonr0(i,2)* gkovr(i,2)
      t23 = t23 + gkonr0(i,2)* gkovr(i,3)

      t31 = t31 + gkonr0(i,3)* gkovr(i,1)
      t32 = t32 + gkonr0(i,3)* gkovr(i,2)
      t33 = t33 + gkonr0(i,3)* gkovr(i,3)
      enddo
      
      !Matrix T
      
      T0mat(1,1) = t11*t11*factor0
      T0mat(1,2) = t11*t21*factor0
      T0mat(1,3) = t11*t31*factor0
      
      T0mat(1,4) = t21*t21*factor0
      T0mat(1,5) = t21*t31*factor0
      T0mat(1,6) = t31*t31*factor0
      
      
      
      T0mat(2,1) = (t11*t12 + t12*t11)*factor0
      T0mat(2,2) = (t11*t22 + t12*t21)*factor0
      T0mat(2,3) = (t11*t32 + t12*t31)*factor0
      
      T0mat(2,4) = (t21*t22 + t22*t21)*factor0
      T0mat(2,5) = (t21*t32 + t22*t31)*factor0
      T0mat(2,6) = (t31*t32 + t32*t31)*factor0
      
      
      T0mat(3,1) = (t11*t13 + t13*t11)*factor0
      T0mat(3,2) = (t11*t23 + t13*t21)*factor0
      T0mat(3,3) = (t11*t33 + t13*t31)*factor0
      
      T0mat(3,4) = (t21*t23 + t23*t21)*factor0
      T0mat(3,5) = (t21*t33 + t23*t31)*factor0
      T0mat(3,6) = (t31*t33 + t33*t31)*factor0
      
      
      T0mat(4,1) = t12*t12*factor0
      T0mat(4,2) = t12*t22*factor0
      T0mat(4,3) = t12*t32*factor0
      
      T0mat(4,4) = t22*t22*factor0
      T0mat(4,5) = t22*t32*factor0
      T0mat(4,6) = t32*t32*factor0
      
      
      
      T0mat(5,1) = (t12*t13 + t13*t12)*factor0
      T0mat(5,2) = (t12*t23 + t13*t22)*factor0
      T0mat(5,3) = (t12*t33 + t13*t32)*factor0
      
      T0mat(5,4) = (t22*t23 + t23*t22)*factor0
      T0mat(5,5) = (t22*t33 + t23*t32)*factor0
      T0mat(5,6) = (t32*t33 + t33*t32)*factor0
      
      
      T0mat(6,1) = t13*t13*factor0
      T0mat(6,2) = t13*t23*factor0
      T0mat(6,3) = t13*t33*factor0
      
      T0mat(6,4) = t23*t23*factor0
      T0mat(6,5) = t23*t33*factor0
      T0mat(6,6) = t33*t33*factor0


      do i=1,6
        do j=1,6
        
        if (abs(T0mat(i,j)).LT.TOL) then
            T0mat(i,j) =  0.d0
        endif  
        
        enddo
      enddo      
              
      
       call kmatrixmultiS10(Weas,6,nEAS,T0mat,6,6,
     * Weasloc,6,nEAS)  

       do i=1,6
        do j=1,nEAS
        if (abs(Weas(i,j)).LE.TOL) then
            Weas(i,j) =  0.d0
        endif  
        
        enddo
      enddo  


      RETURN
      END   
!-------------------------------------------
!---------------------------------------------------------
!---------------------------------------------------------
      subroutine kresidualEAS(FEAS,Weas,stress,nEAS)
      
         !IMPLICIT NONE
  
      INTEGER i,j,k,l,nEAS    
      REAL*8 FEAS(nEAS,1),Weas(6,nEAS), WeasT(nEAS,6),stress(6), stressAux(6,1)

      !initialization
        FEAS(:,1) = 0.d0
      
      do i=1,6
      stressAux(i,1) = 0.d0
      stressAux(i,1) = stress(i) 
      enddo

      call kCompTransposeS10(Weas,nEAS,6,WeasT)

      call kmatrixmultiS10(FEAS,nEAS,1,WeasT,nEAS,6,
     * stressAux,6,1)
     
     
      RETURN
      END   
!---------------------------------------------------------
!---------------------------------------------------------
           
       subroutine kmatricesEAS(KdEAS,KEASEAS,KEASd,Weas,Dmat,Bop,nEAS)

          !IMPLICIT NONE
            
        !       number definition
        
        REAL*8 KdEAS(24,nEAS),KEASEAS(nEAS,nEAS),KEASd(nEAS,24) 
        REAL*8 Weas(6,nEAS), WeasT(nEAS,6), Dmat(6,6), Bop(6,24),BopT(24,6)  
        INTEGER i,j,k,l,nEAS
      !--------------------------  
      !Compute transpose matrices
      !--------------------------  
      call kCompTransposeS10(Weas,nEAS,6,WeasT)

      call kCompTransposeS10(Bop,24,6,BopT)
      !-------------------------- 
      !-------------------------- 
      ! Compute KdEAS
      !--------------------------   
        KdEAS(:,:) = 0.d0     
      call kmatrixtriplS10(KdEAS,24,nEAS,BopT,24,6,
     * Dmat,6,6,Weas,6,nEAS)  
      !-------------------------- 
      ! Compute KEASEAS
      !--------------------------  
       call kmatrixtriplS10(KEASEAS,nEAS,nEAS,WeasT,nEAS,6,
     * Dmat,6,6,Weas,6,nEAS)
      !-------------------------- 
      ! Compute KEASd
      !-------------------------- 
        KEASd(:,:) = 0.d0
       call kmatrixtriplS10(KEASd,nEAS,24,WeasT,nEAS,6,
     * Dmat,6,6,Bop,6,24)
     
                 
            RETURN
            END
!---------------------------------------------------------
!---------------------------------------------------------
       subroutine kEtilde(Etilde,Weas,EASstrain,nEAS)

         !IMPLICIT NONE

        INTEGER i,j,k,l,m,nEAS
        REAL*8 Weas(6,nEAS),EASstrain(nEAS,1),Etilde(6,1) 

        !initialization
        
        Etilde(:,1) = 0.d0
        call  kmatrixmultiS10(Etilde,6,1,Weas,6,nEAS,
     * EASstrain,nEAS,1)

      RETURN
      END
!---------------------------------------------------------
!---------------------------------------------------------
!---------------------------------------------------------
      subroutine kStiffGeomS10Alt(Kgdd,shapef,dshape,stress,ndim,nnode,
     *  ntens)
      
         !IMPLICIT NONE
      
        REAL*8 shapef(nnode)
        REAL*8 dshape(ndim,nnode)
 
        REAL*8 Kgdd(ndim*nnode,ndim*nnode)
        REAL*8 stress(ntens)

        
        INTEGER ndim,nnode,ntens
        INTEGER i,j,k,l,kk,ll
        
        INTEGER node_st, inode
        
	  REAL*8 dNidr,dNids,dNidt
	
	  REAL*8 dNkdr,dNkds,dNkdt 
        REAL*8 Gamma
        
        INTEGER indexcol
        

        Kgdd(:,:) = 0.0d0    
 
        
      node_st=1 !initialization index of rows

      do i=1,nnode
      
            indexcol = 1 !initialization index of columns
            
	  do k=1,nnode

	dNidr =   dshape(1,i)
	dNids =   dshape(2,i)
	dNidt =   dshape(3,i)
	
	dNkdr =   dshape(1,k)
	dNkds =   dshape(2,k) 
	dNkdt =   dshape(3,k)

	Gamma =dNidr*dNkdr*stress(1)+ (dNidr*dNkds+dNids*dNkdr)*stress(2)+
     *    (dNidr*dNkdt + dNidt*dNkdr)*stress(3) + dNids*dNkds*stress(4)+
     *    (dNids*dNkdt + dNidt*dNkds)*stress(5) + dNidt*dNkdt*stress(6)

        Kgdd(node_st,indexcol) =  Gamma
        Kgdd(node_st+1,indexcol+1) =  Gamma
        Kgdd(node_st+2,indexcol+2) =  Gamma
 
                    indexcol = indexcol+3 !update index of columns
                  
                    enddo !end k
                    
	 node_st = node_st+3 !update index of rwos	
	enddo !end i
	
      
        RETURN 
        END
        
!---------------------------------------------------------  
!---------------------------------------------------------
!---------------------------------------------------------
      subroutine kStiffGeomS10AltANS(Kgdd,shapef,dshape,stress,ndim,
     * nnode,ntens,gkovc1q1,gkovc2q2,frq,fsq,dshape1q1,dshape2q2,
     * nANSshear,ANSshearpar,gkovcT,frT,dshapeT,nANStrap,ANStrappar)
          
         !IMPLICIT NONE


        REAL*8 shapef(nnode), transcart(6,6),dshape(ndim,nnode) 
        REAL*8 Kgdd(ndim*nnode,ndim*nnode), stress(ntens),ANSshearpar    
        REAL*8 gkovc1q1(nANSshear,ndim,ndim),gkovc2q2(nANSshear,ndim,ndim) 
        REAL*8 frq(nANSshear),dshape1q1(nANSshear,ndim,nnode),ANStrappar
        REAL*8 fsq(nANSshear),dshape2q2(nANSshear,ndim,nnode)
        REAL*8 gkovcT(nANStrap,ndim,ndim),frT(nANStrap),dshapeT(nANStrap,ndim,nnode)
        INTEGER nANStrap,nANSshear,ndim,nnode,ntens,i,j,k,l,kk,ll
        INTEGER node_st, inode,indexcol  
	  REAL*8 dNidr,dNids,dNidt,dNkdr,dNkds,dNkdt,Gamma, gam_tra(6), gam(6)
        REAL*8 Gamma11,Gamma12,Gamma22,Gamma13,Gamma23,Gamma33
        
        !initialization
        do i=1,nnode*ndim
            do j=1,nnode*ndim
        Kgdd(i,j) = 0.d0    
            enddo
        enddo  
        
        Gamma = 0.d0
        Gamma11 = 0.d0
        Gamma12 = 0.d0
        Gamma13 = 0.d0
        Gamma22 = 0.d0
        Gamma23 = 0.d0
        Gamma33 = 0.d0
      !--------  end of initialization 
      
      node_st=1 !initialization index of rows

      do i=1,nnode
      
            indexcol = 1 !initialization index of columns
            
	  do k=1,nnode

	dNidr =   dshape(1,i)
	dNids =   dshape(2,i)
	dNidt =   dshape(3,i)
	
	dNkdr =   dshape(1,k)
	dNkds =   dshape(2,k) 
	dNkdt =   dshape(3,k)
	

      Gamma11 = dNidr*dNkdr
	Gamma12 = dNidr*dNkds+dNids*dNkdr
	Gamma22 = dNids*dNkds
	
	
	!------- ANS shear
	if (ANSshearpar.ne.0.d0) then
	
	Gamma13 = frq(1)*(dshape1q1(1,1,i)*dshape1q1(1,3,k)) + 
     *          frq(1)*(dshape1q1(1,3,i)*dshape1q1(1,1,k)) + 
     *          frq(2)*(dshape1q1(2,1,i)*dshape1q1(2,3,k)) + 
     *          frq(2)*(dshape1q1(2,3,i)*dshape1q1(2,1,k)) 
     
     	Gamma23 = fsq(1)*(dshape2q2(1,2,i)*dshape2q2(1,3,k)) + 
     *          fsq(1)*(dshape2q2(1,3,i)*dshape2q2(1,2,k)) + 
     *          fsq(2)*(dshape2q2(2,2,i)*dshape2q2(2,3,k)) + 
     *          fsq(2)*(dshape2q2(2,3,i)*dshape2q2(2,2,k)) 
     
     

      else
     
      Gamma13 = dNidr*dNkdt + dNidt*dNkdr
      Gamma23 = dNids*dNkdt + dNidt*dNkds
     
	endif
	
	!------- ANS trap
	if (ANStrappar.ne.0.d0) then
	

	
	Gamma33 = frT(1)*dshapeT(1,3,i)*dshapeT(1,3,k) + 
     *          frT(2)*dshapeT(2,3,i)*dshapeT(2,3,k) +
     *          frT(3)*dshapeT(3,3,i)*dshapeT(3,3,k) +
     *          frT(4)*dshapeT(4,3,i)*dshapeT(4,3,k) 
	
	
	else 
	
	Gamma33 = dNidt*dNkdt
	
	endif
		
		gam_tra(1)=Gamma11
		gam_tra(2)=Gamma12
		gam_tra(3)=Gamma13
		gam_tra(4)=Gamma22
		gam_tra(5)=Gamma23
		gam_tra(6)=Gamma33
		
		!CALL KMATRIXMULTS10TM(gam_tra,6,1,transcart,6,6,
!     *   gam,6,1)
!	Gamma =dNidr*dNkdr*stress(1)+ (dNidr*dNkds+dNids*dNkdr)*stress(2)+
!     *    (dNidr*dNkdt + dNidt*dNkdr)*stress(3) + dNids*dNkds*stress(4)+
!     *    (dNids*dNkdt + dNidt*dNkds)*stress(5) + dNidt*dNkdt*stress(6)
     
     	Gamma =gam_tra(1)*stress(1)+ (gam_tra(2))*stress(2)+
     *    (gam_tra(3))*stress(3) + gam_tra(4)*stress(4)+
     *    (gam_tra(5))*stress(5) + gam_tra(6)*stress(6)

        Kgdd(node_st,indexcol) =  Gamma
        Kgdd(node_st+1,indexcol+1) =  Gamma
        Kgdd(node_st+2,indexcol+2) =  Gamma
 
         indexcol = indexcol+3 !update index of columns
                  
                    enddo !end k
                    
	 node_st = node_st+3 !update index of rwos	
	enddo !end i
	   
      RETURN
      END
     
!---------------------------------------------------------    
!---------------------------------------------------------    
!---------------------------------------------------------    
!---------------------------------------------------
!---------------------------------------------------
      subroutine kCompTransposeS10(A,m,n,AT)  
      
         !IMPLICIT NONE
      
        REAL*8 A(n,m),AT(m,n)            
       INTEGER n,m,i,j,k !indices 
       
        AT(:,:) = 0.d0 !initialization
      
       do i=1,n
        do j=1,m
          AT(j,i) = A(i,j)
          end do
       end do
       
      RETURN
      END
!---------------------------------------------------
c----------------------------------------------
c----------------------------------------------
c----------------------------------------------
! matrix multiplication, Abaqus has this as well directly.----
        subroutine  kmatrixmultiS10(A,a1,a2,B,b1,b2,
     & C,c1,c2)  

c        This subroutine computing 
c        A = B*C
c

        integer a1, a2, b1, b2, c1, c2,i,j,k

        REAL*8 A(a1,a2), B(b1,b2),C(c1,c2), sump    
        !initialization to 0.d0 of the result matrix
        do i=1,a1
            do j=1,a2
        A(i,j) = 0.d0    
            enddo
        enddo   

        if(b2-c1.ne.0.d0) then !It can not be multiplied
        !!write (7,*)'No matrix multiplication'
        end if 

c        operation = factor1*factor2  !operation

        do i=1,b1                
                do j=1,c2
            sump = 0.d0
                        do k=1,c1
                sump=sump+B(i,k)*C(k,j)
                        end do 
                A(i,j)=sump
                end  do 
        end do

        return 
        end 
c---------------------------------        
!---------------------------------Matrix triple multiplications
       subroutine  kmatrixtriplS10(Ares,a1,a2,B,b1,b2,
     & C,c1,c2,Trd,d1,d2)  

c        This subroutine computing 
c        A = B*C*D
c
           !IMPLICIT NONE

        integer a1, a2,b1, b2,c1, c2,d1, d2,i,j,k
        REAL*8 Ares(a1,a2),B(b1,b2), C(c1,c2),Trd(d1,d2),Aux(c1,d2),sump, sump2     
		
        !initialization to 0.d0 of the auxiliar matrix
       do i=1,c1
        do j=1,d2
        Aux(i,j) = 0.d0
        enddo
       enddo
       
       do i=1,a1
        do j=1,a2
        Ares(i,j) = 0.d0
         enddo
       enddo

       if(b2-c1.ne.0.d0) then !It can not be multiplied
           !!!write(7,*)'No matrix multiplication'
       end if 
       
       if(c2-d1.ne.0.d0) then !It can not be multiplied
           !!!write(7,*)'No matrix multiplication'
       end if 
     
      call kmatrixmultiS10(Aux,c1,d2,C,c1,c2,Trd,d1,d2)
        
      call kmatrixmultiS10(Ares,b1,d2,B,b1,c1,Aux,c1,d2)
       

        return 
        end 
!---------------------------------------------------------
!---------------------------------------------------------
!---------------------------------------------------------
  
	subroutine InverMatrixS10(A,Ainv,N) 
	
	   !IMPLICIT NONE

	INTEGER N,indx(N),i, j, m1, m2
	REAL*8 A(N,N), Ainv(N,N),y(n,n),col(N) !original matrix
	REAL*8 :: code = 0.d0
	REAL*8 ::  d = 1.d0
	
	!initialization of the inverse
	Ainv(:,:) = 0.0d0
	do i=1,n
		do j=1,n
	Ainv(i,j) = A(i,j) 
		end do
	end do

	!INITIALIZATION OF INDX
	indx(:) = 0.d0
	
	y(:,:)=0.d0
	do i=1,n !Set up identity matrix. 
		y(i,i)=1.d0
	end do 
	
	call ludcmpS10(Ainv,n,indx,d,code) !Decompose the matrix just once.

			do  j=1,n !Find inverse by columns.
				call lubksbS10(Ainv,n,indx,y(1,j))
			end do 
		Ainv(:,:) = y(:,:)
	return 
	end
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
	Subroutine LUDCMPS10(A,N,INDX,D,CODE)
	
	REAL*8,PARAMETER :: NMAX=100,TINY=1.5D-16
	REAL*8  AMAX,DUM, sum2,CODE, A(N,N),VV(NMAX)
	INTEGER N, INDX(N),I, J, K, IMAX
	

	D=1; CODE=0

	DO I=1,N
		 AMAX=0.d0
			DO J=1,N
			IF (DABS(A(I,J)).GT.AMAX) AMAX=DABS(A(I,J))
			END DO ! j loop
		IF(AMAX.LT.TINY) THEN
			CODE = 1
				RETURN
		END IF
		  VV(I) = 1.d0 / AMAX
	END DO ! i loop

	DO J=1,N
			DO I=1,J-1
			sum2 = A(I,J)
				DO K=1,I-1
				sum2 = sum2 - A(I,K)*A(K,J) 
				END DO ! k loop

			A(I,J) = sum2

			END DO ! i loop

		AMAX = 0.d0

		DO I=J,N
			sum2 = A(I,J)
				DO K=1,J-1
				sum2 = sum2 - A(I,K)*A(K,J) 
				END DO ! k loop
			A(I,J) = sum2
		DUM = VV(I)*DABS(sum2)
			IF(DUM.GE.AMAX) THEN
			IMAX = I
			AMAX = DUM
			END IF
		END DO ! i loop  
   
		IF(J.NE.IMAX) THEN
      
		DO K=1,N
				DUM = A(IMAX,K)
				A(IMAX,K) = A(J,K)
				A(J,K) = DUM
		END DO ! k loop
			
		
			D = -D
		
		VV(IMAX) = VV(J)
		END IF

			
			INDX(J) = IMAX
      
			IF(DABS(A(J,J)).LT. TINY) A(J,J) = TINY

		IF(J.NE.N) THEN
			DUM = 1.d0 / A(J,J)
			DO I=J+1,N
				A(I,J) = A(I,J)*DUM
			END DO ! i loop
		END IF 
      END DO ! j loop

      RETURN
      END
!------------------------------------------------------------------
!-------------------------------------------------------------------
	Subroutine LUBKSBS10(A,N,INDX,B)


	INTEGER N,INDX(N), II, I, J, LL
	REAL*8 A(N,N),B(N),sum2 

		II = 0

		DO I=1,N
			LL = INDX(I)
			sum2 = B(LL)
			B(LL) = B(I)
		IF(II.NE.0) THEN
			 DO J=II,I-1
			sum2 = sum2 - A(I,J)*B(J)
          END DO ! j loop
			 ELSE IF(sum2.NE.0.d0) THEN
				 II = I
					END IF
			B(I) = sum2
		END DO ! i loop

		 DO I=N,1,-1
			sum2 = B(I)
		IF(I.LT.N) THEN
       DO J=I+1,N
			sum2 = sum2 - A(I,J)*B(J)
         END DO ! j loop
      END IF
       B(I) = sum2 / A(I,I)
       END DO ! i loop

        RETURN 
        END
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!--------------------------------------------------------- 
      subroutine kshapefANSshear(r,s,frq,fsq)

      REAL*8 frq(2),fsq(2),r,s
      !initialization
      frq(:) = 0.d0
      fsq(:) = 0.d0

      frq(1) = 0.5d0 *(1.d0 -s)
      frq(2) = 0.5d0 *(1.d0 +s)

      fsq(1) = 0.5d0 *(1.d0 -r)
      fsq(2) = 0.5d0 *(1.d0 +r)
      RETURN
      END
!---------------------------------------------------------
!---------------------------------------------------------       

      subroutine kANScolpointsShear(gkovr1q1,gkovc1q1,gkonr1q1,gkonc1q1,
     * gmkovr1q1,gmkovc1q1,gmkonr1q1,gmkonc1q1,detr1,detc1,shapef1q1,
     * dshape1q1,gkovr2q2,gkovc2q2,gkonr2q2,gkonc2q2,gmkovr2q2,
     * gmkovc2q2,gmkonr2q2,gmkonc2q2,detr2,detc2,shapef2q2,dshape2q2,
     * xref,xcur,ndim,nnode,nANSshear)
     
         !IMPLICIT NONE


        
        INTEGER i,j,k,l,m,kk,pp,qq,mm,iq
        
      REAL*8 xref(ndim,nnode) !REFERENCE CONFIG COORDINATES
      REAL*8 xcur(ndim,nnode) !CURRENT CONFIG COORDINATES
        
        !ANS parameter 
        INTEGER nANSshear,ndim,nnode
        
        !ANS variables
        REAL*8 gkovr1q1(nANSshear,ndim,ndim),gkovc1q1(nANSshear,ndim,ndim)
        REAL*8 gkonr1q1(nANSshear,ndim,ndim),gkonc1q1(nANSshear,ndim,ndim)
        REAL*8 gmkovr1q1(nANSshear,ndim,ndim),gmkovc1q1(nANSshear,ndim,ndim)
        REAL*8 gmkonr1q1(nANSshear,ndim,ndim),gmkonc1q1(nANSshear,ndim,ndim)
        REAL*8 detr1(nANSshear), detc1(nANSshear),detr2(nANSshear),detc2(nANSshear)
        REAL*8 shapef1q1(nANSshear,nnode), dshape1q1(nANSshear,ndim,nnode)
        REAL*8 gkovr2q2(nANSshear,ndim,ndim),gkovc2q2(nANSshear,ndim,ndim)
        REAL*8 gkonr2q2(nANSshear,ndim,ndim),gkonc2q2(nANSshear,ndim,ndim)
        REAL*8 gmkovr2q2(nANSshear,ndim,ndim),gmkovc2q2(nANSshear,ndim,ndim)
        REAL*8 gmkonr2q2(nANSshear,ndim,ndim),gmkonc2q2(nANSshear,ndim,ndim)
        REAL*8 shapef2q2(nANSshear,nnode),dshape2q2(nANSshear,ndim,nnode)
        REAL*8 shapef1qi(nnode),dshape1qi(ndim,nnode),gkovr1qi(ndim,ndim)
        REAL*8 gkonr1qi(ndim,ndim),gmkovr1qi(ndim,ndim),gmkonr1qi(ndim,ndim)
        REAL*8 gkovc1qi(ndim,ndim),detr1i,gkonc1qi(ndim,ndim),gmkovc1qi(ndim,ndim)
        REAL*8 gmkonc1qi(ndim,ndim),detc1i,shapef2qi(nnode),dshape2qi(ndim,nnode)
        REAL*8 gkovr2qi(ndim,ndim),gkonr2qi(ndim,ndim), gmkovr2qi(ndim,ndim)
        REAL*8 gmkonr2qi(ndim,ndim),detr2i, gkovc2qi(ndim,ndim), gkonc2qi(ndim,ndim)
        REAL*8 gmkovc2qi(ndim,ndim),gmkonc2qi(ndim,ndim),detc2i
       REAL*8 xr1(2),xs1(2), xr2(2),xs2(2) 
       
       !initialization of variables

       xr1(:) = 0.d0
       xs1(:) = 0.d0
       
       xr2(:) = 0.d0
       xs2(:) = 0.d0

            detr1(:) =  0.d0
            detc1(:) =  0.d0
            
            detr2(:) =  0.d0
            detc2(:) =  0.d0
       
         gkovr1q1(:,:,:) = 0.d0
         gkovc1q1(:,:,:) = 0.d0
         gkonr1q1(:,:,:) = 0.d0
         gkonc1q1(:,:,:) = 0.d0
        
         gmkovr1q1(:,:,:) = 0.d0
         gmkovc1q1(:,:,:) = 0.d0
         gmkonr1q1(:,:,:) = 0.d0
         gmkonc1q1(:,:,:) = 0.d0
         
         gkovr2q2(:,:,:) = 0.d0
         gkovc2q2(:,:,:) = 0.d0
         gkonr2q2(:,:,:) = 0.d0
         gkonc2q2(:,:,:) = 0.d0
      
         gmkovr2q2(:,:,:) = 0.d0
         gmkovc2q2(:,:,:) = 0.d0
         gmkonr2q2(:,:,:) = 0.d0
         gmkonc2q2(:,:,:) = 0.d0
       shapef1q1(:,:) =  0.d0
       shapef2q2(:,:) =  0.d0
        dshape1q1(:,:,:) =  0.d0
        dshape2q2(:,:,:) =  0.d0 
       !--------------------------
        xr1(1) =  0.d0
        xs1(1) = -1.d0
        
        xr1(2) = 0.d0
        xs1(2) = 1.d0
        
        xr2(1) = -1.d0 
        xs2(1) =  0.d0
        
        xr2(2) =  1.d0   
        xs2(2) =  0.d0  
      !--------------------------
      
       do iq =1,nANSshear
       !Compute shape functions and derivatives
      call kShapeFunctS10(xr1(iq),xs1(iq),0.d0,ndim,nnode,
     * shapef1qi,dshape1qi)
     
            do  i=1,nnode
            shapef1q1(iq,i) = shapef1qi(i)
            enddo
            
            do i=1,ndim
                do j=1,nnode
            dshape1q1(iq,i,j) =  dshape1qi(i,j) 
                enddo
            enddo
          
      !Compute metrics at the reference and current configurations at collocation points
       call ks8metShellBody10(xref,shapef1qi,dshape1qi,0.d0,nnode,ndim, 
     * gkovr1qi,gkonr1qi,gmkovr1qi,gmkonr1qi,detr1i)
     
       call ks8metShellBody10(xcur,shapef1qi,dshape1qi,0.d0,nnode,ndim, 
     * gkovc1qi,gkonc1qi,gmkovc1qi,gmkonc1qi,detc1i)


            
         do i=1,ndim
           do j=1,ndim
            gkovr1q1(iq,i,j) =   gkovr1qi(i,j)   
            gkovc1q1(iq,i,j) =   gkovc1qi(i,j) 
    
            gkonr1q1(iq,i,j) =   gkonr1qi(i,j)
            gkonc1q1(iq,i,j) =   gkonc1qi(i,j)
    
            gmkovr1q1(iq,i,j) = gmkovr1qi(i,j)
            gmkovc1q1(iq,i,j) = gmkovc1qi(i,j)
    
            gmkonr1q1(iq,i,j) = gmkonr1qi(i,j)
            gmkonc1q1(iq,i,j) = gmkonc1qi(i,j)
    
            detr1(iq) = detr1i
            detc1(iq) = detc1i
            enddo
         enddo
            
!--------------------------------------------------------------------            
                     
       !Compute shape functions and derivatives
      call kShapeFunctS10(xr2(iq),xs2(iq),0.d0,ndim,nnode,
     * shapef2qi,dshape2qi)
     
            do  i=1,nnode
            shapef2q2(iq,i) = shapef2qi(i)
            enddo
            
            do i=1,ndim
                do j=1,nnode
            dshape2q2(iq,i,j) =  dshape2qi(i,j) 
                enddo
            enddo
                        
      !Compute metrics at the reference and current configurations at collocation points
       call ks8metShellBody10(xref,shapef2qi,dshape2qi,0.d0,nnode,ndim, 
     * gkovr2qi,gkonr2qi,gmkovr2qi,gmkonr2qi,detr2i)
     
       call ks8metShellBody10(xcur,shapef2qi,dshape2qi,0.d0,nnode,ndim, 
     * gkovc2qi,gkonc2qi,gmkovc2qi,gmkonc2qi,detc2i)
              
         do i=1,ndim
           do j=1,ndim
            gkovr2q2(iq,i,j) =   gkovr2qi(i,j)   
            gkovc2q2(iq,i,j) =   gkovc2qi(i,j) 
    
            gkonr2q2(iq,i,j) =   gkonr2qi(i,j)
            gkonc2q2(iq,i,j) =   gkonc2qi(i,j)
    
            gmkovr2q2(iq,i,j) = gmkovr2qi(i,j)
            gmkovc2q2(iq,i,j) = gmkovc2qi(i,j)
    
            gmkonr2q2(iq,i,j) = gmkonr2qi(i,j)
            gmkonc1q1(iq,i,j) = gmkonc2qi(i,j)
    
            detr2(iq) = detr2i
            detc2(iq) = detc2i
            enddo
         enddo               
               
      enddo !end loop over the collocation points
      
             
      RETURN
      END  
!--------------------------------------------------------- 
!--------------------------------------------------------- 
      subroutine kANSshearS10(Bop,ndim,nnode,ntens,ndof,gkovc1q1,
     * gkovc2q2,frq,fsq,dshape1q1,dshape2q2,nANSshear)
     
         !IMPLICIT NONE
        
        INTEGER ndim,nnode,ntens,ndof, nANSshear, i, j, node_st,innode
        
      REAL*8 Bop(ntens,ndim*nnode), Bop_curv(ntens,ndim*nnode)
      REAL*8 BopAux(ntens,ndim*nnode),gkovc1q1(nANSshear,ndim,ndim)   
      REAL*8 gkovc2q2(nANSshear,ndim,ndim), frq(nANSshear), fsq(nANSshear)
      REAL*8 dshape1q1(nANSshear,ndim,nnode), dshape2q2(nANSshear,ndim,nnode)
      
      !initialization

      BopAux(:,:) = 0.d0
	  BopAux(:,:) =  Bop(:,:)
 
      BopAux(3,:) = 0.d0
      BopAux(5,:) = 0.d0
      
        
       node_st = 1
      do innode = 1,nnode
      BopAux (3,node_st) = BopAux(3,node_st)+ 
     * frq(1)*dshape1q1(1,1,innode)*gkovc1q1(1,1,3) + 
     * frq(1)*dshape1q1(1,3,innode)*gkovc1q1(1,1,1) + 
     * frq(2)*dshape1q1(2,1,innode)*gkovc1q1(2,1,3) + 
     * frq(2)*dshape1q1(2,3,innode)*gkovc1q1(2,1,1) 

      
       BopAux (3,node_st+1) = BopAux(3,node_st+1)+ 
     * frq(1)*dshape1q1(1,1,innode)*gkovc1q1(1,2,3) + 
     * frq(1)*dshape1q1(1,3,innode)*gkovc1q1(1,2,1) + 
     * frq(2)*dshape1q1(2,1,innode)*gkovc1q1(2,2,3) + 
     * frq(2)*dshape1q1(2,3,innode)*gkovc1q1(2,2,1)  

      
       BopAux (3,node_st+2) = BopAux(3,node_st+2)+ 
     * frq(1)*dshape1q1(1,1,innode)*gkovc1q1(1,3,3) + 
     * frq(1)*dshape1q1(1,3,innode)*gkovc1q1(1,3,1) + 
     * frq(2)*dshape1q1(2,1,innode)*gkovc1q1(2,3,3) + 
     * frq(2)*dshape1q1(2,3,innode)*gkovc1q1(2,3,1) 
 
     
      BopAux (5,node_st) = BopAux(5,node_st) +
     * fsq(1)*dshape2q2(1,2,innode)*gkovc2q2(1,1,3) + 
     * fsq(1)*dshape2q2(1,3,innode)*gkovc2q2(1,1,2) +
     * fsq(2)*dshape2q2(2,2,innode)*gkovc2q2(2,1,3) + 
     * fsq(2)*dshape2q2(2,3,innode)*gkovc2q2(2,1,2) 
     
     
       BopAux (5,node_st+1) = BopAux(5,node_st+1) +
     * fsq(1)*dshape2q2(1,2,innode)*gkovc2q2(1,2,3) + 
     * fsq(1)*dshape2q2(1,3,innode)*gkovc2q2(1,2,2) +
     * fsq(2)*dshape2q2(2,2,innode)*gkovc2q2(2,2,3) + 
     * fsq(2)*dshape2q2(2,3,innode)*gkovc2q2(2,2,2) 

     
       BopAux (5,node_st+2) = BopAux(5,node_st+2) +
     * fsq(1)*dshape2q2(1,2,innode)*gkovc2q2(1,3,3) + 
     * fsq(1)*dshape2q2(1,3,innode)*gkovc2q2(1,3,2) +
     * fsq(2)*dshape2q2(2,2,innode)*gkovc2q2(2,3,3) + 
     * fsq(2)*dshape2q2(2,3,innode)*gkovc2q2(2,3,2)   

      node_st  = node_st  + ndim
      enddo
 
		Bop(:,:) =   BopAux(:,:)
      
      RETURN
      END
!--------------------------------------------------------- 
!---------------------------------------------------------       
!---------------------------------------------------------       
!---------------------------------------------------------      
      subroutine kshapefANStrap(r,s,frT)   

      REAL*8 frT(4), r,s
      frT(:) =  0.d0
      frT(1) = ((1.d0+r)*(1.d0+s))/4.d0
	frT(2) = ((1.d0-r)*(1.d0+s))/4.d0
	frT(3) = ((1.d0-r)*(1.d0-s))/4.d0
	frT(4) = ((1.d0+r)*(1.d0-s))/4.d0
      
      RETURN
      END
!--------------------------------------------------------- 
!--------------------------------------------------------- 
      subroutine kANStrapS10(Bop,ndim,nnode,ntens,ndof,gkovcT,
     * frT,dshapeT,nANStrap)
     
              !IMPLICIT NONE
       
        
        INTEGER i,j,node_st,innode,jq
      INTEGER ndim,nnode,ntens,ndof,nANStrap
        
      REAL*8 Bop(6,ndim*nnode),Bop_curv(6,ndim*nnode),BopAux(6,ndim*nnode) 
      REAL*8 frT(nANStrap), gkovcT(nANStrap,ndim,ndim),dshapeT(nANStrap,ndim,nnode) 
            !initialization

      BopAux(:,:) = 0.d0


		BopAux(:,:) =  Bop(:,:) 
       node_st = 1
      do innode = 1,nnode
        do jq=1,nANStrap
		
      BopAux(6,node_st) =   BopAux(6,node_st)   + 
     *   frT(jq)*dshapeT(jq,3,innode)*gkovcT(jq,1,3)
     
      BopAux(6,node_st+1) = BopAux(6,node_st+1) +
     *   frT(jq)*dshapeT(jq,3,innode)*gkovcT(jq,2,3)
      
      BopAux(6,node_st+2) = BopAux(6,node_st+2) +
     *  frT(jq)*dshapeT(jq,3,innode)*gkovcT(jq,3,3)
        enddo
      
       node_st  =  node_st  + ndim
      enddo
	  
      Bop(:,:) =   BopAux(:,:) 
      
      RETURN
      END
!---------------------------------------------------------    

      subroutine kMatLawS10(Cmat,gmkovc,gmkovr,gmkonr,gkonr,gkovr,
     * DefGrad,DetF,DefGradInv, emod,enu,stress,strain, 
     *	 ntens,kFreeEner,kFreeEner_pos,kFreeEner_neg)
        
        IMPLICIT NONE
      INTEGER ntens, Nota(3,3), i, j, k, kk,mm, pp, qq,l

      REAL*8 Cmat(ntens,ntens), transcart(6,6),EGL_cart(3,3),EGL_princ(3)
      REAL*8 stress(ntens), Cmat_Cart(6,6), emod,enu,kFreeEner,Edeviotoric(3)
      REAL*8 TraceE,Norm,  TraceE_pos, TraceE_neg,kFreeEner_pos,kFreeEner_neg
      REAL*8 strain(ntens), gmkovc(3,3), gmkovr(3,3),gmkonr(3,3)
      REAL*8 Ev(ntens), El(ntens), EGL(3,3), C(3,3,3,3), gkonr(3,3)
      REAL*8 Strain_cart(ntens), Strain_norCart(ntens)
      REAL*8 DefGrad(3,3), DetF, DefGradInv(3,3),DefGradInvT(3,3)
      REAL*8  c1, c2,bk,eg,elam,sum2,gkovr(3,3)
	   REAL*8 Euler_cart(3,3),DefGradT(3,3),cauchy(3,3), cauchy1(3,3)
	  REAL*8 StressL(3,3), EGL_cart1(3,3)
      REAL*8 sig_pos(3), tr_sig, tr_sig_pos,ene_crk,stress_princ(3), stress_cart(3,3)


		

        Data Nota/1,2,3,             
     *            2,4,5,
     *            3,5,6/  
      
      !initialization
      EGL(:,:) = 0.d0
      C(:,:,:,:) = 0.d0
      stress(:) = 0.d0
      strain(:) = 0.d0
      Ev(:) = 0.d0
      Cmat(:,:) = 0.d0
	  EGL_cart1(:,:)=0.0d0
	  !EGL(:,:)=0.0d0
	  Euler_cart(:,:)=0.0d0
	  EGL_princ(:)=0.0d0
	  stressL(:,:)=0.0d0
	  stress_cart(:,:)=0.0d0
	  cauchy1(:,:)=0.0d0
	  cauchy(:,:)=0.0d0
	  stress_princ(:)=0.0d0
	  

c     Green-Lagrange strain tensor  in vector form
      Ev(1) = (0.5d0) * (gmkovc(1,1) - gmkovr(1,1) )  !E11
      Ev(2) = (0.5d0) * (gmkovc(1,2) - gmkovr(1,2) )  !E12
      Ev(3) = (0.5d0) * (gmkovc(1,3) - gmkovr(1,3) )  !E13
      Ev(4) = (0.5d0) * (gmkovc(2,2) - gmkovr(2,2) )  !E22
      Ev(5) = (0.5d0) * (gmkovc(2,3) - gmkovr(2,3) )  !E23
      Ev(6) = (0.5d0) * (gmkovc(3,3) - gmkovr(3,3) )  !E33
      
	  strain(:)=Ev(:)

      EGL(1,1)  = Ev(1)
      EGL(1,2)  = Ev(2) 
      EGL(1,3)  = Ev(3)
      
      EGL(2,1) = Ev(2)
      EGL(2,2) = Ev(4)
      EGL(2,3) = Ev(5)
      
      EGL(3,1) = Ev(3)
      EGL(3,2) = Ev(5)
      EGL(3,3) = Ev(6)
	  
	  
		
! MATERIAL KIRCHHOFF-SAINT-Venant
      c1 = (emod*enu)/((1.0d0+enu)*(1.0d0-2.0d0*enu))
	 c2 = (emod)/(2.0d0*(1.0d0+enu))

    		     
	   bk=emod/3.0d0/(1.0d0-2.0d0*enu)
	    eg=emod/2.0d0/(1.0d0+enu)
	     elam=bk-2.0d0*eg/3.0d0
		 
      do i=1,3
        do  j=1,3
            do  k=1,3
                do  l=1,3 
      C(i,j,k,l) = c1*gmkonr(i,j)*gmkonr(k,l) +
     *   c2*(gmkonr(i,k)*gmkonr(j,l) + gmkonr(i,l)*gmkonr(k,j))
                enddo
            enddo
        enddo
      enddo
      
        do kk=1,3
           do mm=1,3
             do pp=1,3
               do qq=1,3
               Cmat(nota(kk,mm),nota(pp,qq))=C(kk,mm,pp,qq)
               end do
             end do
           end do
         end do


      !normal strains
      El(1) =Ev(1);
      El(4) =Ev(4);
      El(6) =Ev(6);
      !shear strain 
      El(2) =Ev(2)*2.d0
      El(3) =Ev(3)*2.d0
      El(5) =Ev(5)*2.d0
	  
 
	 ! Since CMAT and strains are already in Cartersian. the stress is also in caretsian 
      do i=1,ntens
        sum2 =0.d0
        do k=1,ntens
            sum2 = sum2 + Cmat(k,i)*El(k) 
        enddo
        stress(i) = sum2  ! This is in curvilinear as well  
      enddo
      kFreeEner = 0.0d0
	  
      do i=1,6
        kFreeEner = kFreeEner + 0.5d0*stress(i)*El(i)
      enddo


		
		
		stressL(1,1)=stress(1)
		stressL(1,2)=stress(4)
		stressL(1,3)=stress(5)
		
		stressL(2,1)=stress(4)
		stressL(2,2)=stress(2)
		stressL(2,3)=stress(6)
		
		stressL(3,1)=stress(5)
		stressL(3,2)=stress(6)
		stressL(3,3)=stress(3)
		

		
		!convert the stress into cartisean
		call s8_kon_cuca(stressL,stress_cart,gkovr) 
		!compute the PK2 in Cartesian frame 
		
		      !compute the transpose of the deformation gradient
		      call kCompTransposeS10(DefGrad,3,3,DefGradT)

		
		!------------------------------stress------------------------------		
				      call kmatrixmultiS10(cauchy1,3,3,DefGrad,
     * 3,3,stress_cart,3,3)
	 
	 		      call kmatrixmultiS10(cauchy,3,3,cauchy1,
     * 3,3,DefGradT,3,3)
	 
!            You can  compute the Cauchy stress as follows	 

			cauchy(:,:)= cauchy(:,:)/DetF
	

		
		!-----------------------For stresses
!for streses
			CALL SPRINC(cauchy,stress_princ,2,3,3) ! Compute principle stress
		
		sig_pos= 0.5d0*(abs(stress_princ)+stress_princ)
			
			
				
		tr_sig=stress_princ(1)+stress_princ(2)+stress_princ(3)
		tr_sig_pos=0.5d0*(abs(tr_sig)+tr_sig)   
                
		
				! Energy driving the crack growth
		ene_crk = 0.5d0*(1+enu)*(sig_pos(1)**2+sig_pos(2)**2+sig_pos(3)**2)/emod -
     * 0.5d0*enu*(tr_sig_pos**2)/emod 
		!write(7,*)'ene_crk', ene_crk
		kFreeEner_pos=ene_crk
		kFreeEner_neg=0.0d0	
		

		
      RETURN
      END

!---------------------------------------------------------
!------------------------------------------------------------------------	  
       subroutine s8_kov_cuca(EGL,Strain_cart,gkonc)
	  
	  implicit none
!       number definition
  

        INTEGER i,j
        REAL*8 EGL(3,3), Tcur(3,3),Strain_cart(3,3),gkonc(3,3)
		REAL*8 EGL_cart(6)
		
	
		    Tcur(:,:)=0.d0

	
		    Tcur(:,:)=EGL(:,:)


		    Strain_cart(:,:)=0.0d0
	

		do i=1,3
		  do j=1,3
		    Strain_cart(1,1)=Strain_cart(1,1)+gkonc(1,i)*gkonc(1,j)*Tcur(i,j)
		    Strain_cart(1,2)=Strain_cart(1,2)+gkonc(1,i)*gkonc(2,j)*Tcur(i,j)
		    Strain_cart(1,3)=Strain_cart(1,3)+gkonc(1,i)*gkonc(3,j)*Tcur(i,j)
		    Strain_cart(2,2)=Strain_cart(2,2)+gkonc(2,i)*gkonc(2,j)*Tcur(i,j)
		    Strain_cart(2,3)=Strain_cart(2,3)+gkonc(2,i)*gkonc(3,j)*Tcur(i,j)
		    Strain_cart(3,3)=Strain_cart(3,3)+gkonc(3,i)*gkonc(3,j)*Tcur(i,j)			
		  end do
		end do
		
		!Strain_cart(1,2)=Strain_cart(1,2)*2.0d0
		!!Strain_cart(1,3)=Strain_cart(1,3)*2.0d0
		!Strain_cart(3,1)=Strain_cart(3,1)*2.0d0
		!Strain_cart(2,1)=Strain_cart(2,1)*2.0d0
		!Strain_cart(2,3)=Strain_cart(2,3)*2.0d0
		!Strain_cart(3,2)=Strain_cart(3,2)*2.0d0
		
		!EGL_cart(1)=Strain_cart(1,1)
		!EGL_cart(4)=Strain_cart(1,2)*2.0d0
		!EGL_cart(5)=Strain_cart(1,3)*2.0d0
		!EGL_cart(2)=Strain_cart(2,2)
		!EGL_cart(6)=Strain_cart(2,3)*2.0d0
		!EGL_cart(3)=Strain_cart(3,3)
		
		return
		end
!---------------------------------------------------------
!------------------------------------------------------------------------	  
       subroutine s8_kon_cuca(PK2curv,PK2_cart,gkovr)
	  
	  implicit none
!       number definition
  

        INTEGER i,j
        REAL*8 PK2curv(3,3), Tcur(3,3),PK2_cart(3,3),gkovr(3,3)
		 
		
	
		    Tcur(:,:)=0.d0
		    PK2_cart(:,:)=0.0d0

		    Tcur(:,:)=PK2curv(:,:)


		do i=1,3
		  do j=1,3
		    PK2_cart(1,1)=PK2_cart(1,1)+gkovr(1,i)*gkovr(1,j)*Tcur(i,j)
		    PK2_cart(1,2)=PK2_cart(1,2)+gkovr(1,i)*gkovr(2,j)*Tcur(i,j)
		    PK2_cart(1,3)=PK2_cart(1,3)+gkovr(1,i)*gkovr(3,j)*Tcur(i,j)
		    PK2_cart(2,2)=PK2_cart(2,2)+gkovr(2,i)*gkovr(2,j)*Tcur(i,j)
		    PK2_cart(2,3)=PK2_cart(2,3)+gkovr(2,i)*gkovr(3,j)*Tcur(i,j)
		    PK2_cart(3,3)=PK2_cart(3,3)+gkovr(3,i)*gkovr(3,j)*Tcur(i,j)			
		  end do
		end do
		
	
		
		return
		end
!---------------------------------------------------------




       subroutine kStiffMatrixS10(Kdd,Fd,nnode,ndim,ndof,dvol,Dmat,
     *  stress,Bop,ntens)

         !IMPLICIT NONE


        INTEGER nnode,ndim,ndof,ntens,iforce,jforce,istiff,jstiff,i,j,k,l,m,n
        REAL*8 Kdd(ndim*nnode,ndim*nnode),Fd(ndim*nnode),Bop(ntens,ndim*nnode)
        REAL*8 Dmat(ntens,ntens), stress(ntens), db(ntens), sum,dvol, sumforce

		Fd(:)=0.d0
			Kdd(:,:)=0.d0

      do j=1,ndim*nnode
        do k=1,ntens
            db(k) =  0.d0
                do l=1,ntens
                db(k) = db(k) + Dmat(k,l)*Bop(l,j)*dvol
                enddo
        enddo

        do i=1,ndim*nnode
            sum = 0.d0
            do m=1,ntens
            sum = sum + Bop(m,i)*db(m)
            enddo
            Kdd(i,j) = Kdd(i,j) + sum
        enddo
      enddo

      do jforce =1,nnode*ndim
                  sumforce= 0.d0
                        do iforce = 1,ntens
       sumforce = sumforce + Bop(iforce,jforce)*stress(iforce)*dvol
                        end do
                 Fd(jforce) = sumforce
      end do

      RETURN
      END
!---------------------------------------------------------
      subroutine kGradShapeFunctionS10(kGradPhi,Boppf,
     * dshape,gkonr,nnode,ndim)

         !IMPLICIT NONE


        REAL*8 kGradPhi(ndim,nnode), Boppf(ndim,nnode),dshape(ndim,nnode), gkonr(ndim,ndim)
        INTEGER i,j,k,l,m,ndim,nnode
        !initialization

        kGradPhi(:,:) =  0.d0
        Boppf(:,:) = 0.d0
 

      do i=1,nnode
      kGradPhi(1,i) = dshape(1,i)*gkonr(1,1) +
     *                  dshape(2,i)*gkonr(1,2) + dshape(3,i)*gkonr(1,3)

      kGradPhi(2,i) = dshape(1,i)*gkonr(2,1) +
     *                  dshape(2,i)*gkonr(2,2) + dshape(3,i)*gkonr(2,3)

      kGradPhi(3,i) = dshape(1,i)*gkonr(3,1) +
     *                  dshape(2,i)*gkonr(3,2) + dshape(3,i)*gkonr(3,3)
      enddo

      do  i=1,ndim
        do j=1,nnode
        Boppf(i,j) =kGradPhi(i,j)
        enddo
      enddo

      RETURN
      END
!---------------------------------------------------------
!---------------------------------------------------------

      subroutine kPhaseInterpS10(kPhinterp,kGradPhi,shapef,Boppf,
     * ndim,nnode,pf_v,kFreeEner_pos,kFreeEner_neg)

          !IMPLICIT NONE

        REAL*8 kPhinterp,kGradPhi(ndim), pf_v(nnode), shapef(nnode), Boppf(ndim,nnode)
        INTEGER ndim,nnode, i,j,k,l,m,n !index

        !initialization
        kPhinterp =  0.d0
        kGradPhi(:) =  0.d0
   

        do i=1,nnode
        kPhinterp = kPhinterp + shapef(i)*pf_v(i)
        enddo
		
		if (kPhinterp .gt. 1.d0) then
		kPhinterp=1.d0
		endif

       do i=1,ndim
        do j=1,nnode
            kGradPhi(i) = kGradPhi(i) + Boppf(i,j)*pf_v(j)
        enddo
       enddo

      RETURN
      END

!---------------------------------------------------------
!---------------------------------------------------------

      subroutine kPhaseFieldResS10(ppf,kPhinterp,kGradPhi,KFreeEner,Gc,Ls0,
     *    shapef,Boppf,nnode,ndim,factor,Fat_deg)

          !IMPLICIT NONE


       INTEGER i,j,k,l,nnode,ndim
       REAL*8 ppf(nnode,1),fs1(nnode),fs2(nnode),fs3(nnode),Fat_deg
       REAL*8 kPhinterp,kGradPhi(ndim,1),KFreeEner,Gc,Ls0,factor
       REAL*8 shapef(nnode), Boppf(ndim,nnode)


       !Initialization
       
       ppf(:,1) = 0.d0
       fs1(:) = 0.d0
       fs2(:) = 0.d0
       fs3(:) = 0.d0
  

       do i=1,nnode
       fs1(i) =Fat_deg*Gc*Ls0*(Boppf(1,i)*kGradPhi(1,1)+
     *                Boppf(2,i)*kGradPhi(2,1)+
     *                Boppf(3,i)*kGradPhi(3,1))
	   fs2(i) = (Fat_deg*Gc/Ls0 + 2.d0*KFreeEner)*kPhinterp*shapef(i)
	   fs3(i)= 2.d0*shapef(i)*KFreeEner
	   ppf(i,1) = (fs1(i) + fs2(i) - fs3(i))*factor
       enddo

      RETURN
      END
!---------------------------------------------------------
!---------------------------------------------------------
!---------------------------------------------------------
      subroutine  kPhaseFieldStiffS10(spf,KFreeEner,Gc,Ls0,
     *    shapef,Boppf,nnode,ndim,factor,Fat_deg)

         !IMPLICIT NONE


       INTEGER i,j,k,l,ik,nnode,ndim
       REAL*8 spf(nnode,nnode),kss1(nnode,nnode),KFreeEner,Gc, Ls0
       REAL*8 shapef(nnode),Fat_deg, Boppf(ndim,nnode), BoppfT(nnode,ndim)
       REAL*8 kss2(nnode,nnode),factor

 
       !initialization
       spf(:,:) = 0.d0
       kss1(:,:) = 0.d0
       kss2(:,:) = 0.d0     

      !transpose Bmat PS
      call kCompTransposeS10(Boppf,nnode,ndim,BoppfT)
      !kCompTranspose(A,m,n,AT)  A(n,m) , At(m,n)

      call kmatrixmultiS10(kss2,nnode,nnode,BoppfT,
     * nnode,ndim,Boppf,ndim,nnode)
      ! kmatrixmulti(A,a1,a2,B,b1,b2,C,c1,c2)  A = B*C

      do i=1,nnode
        do j=1,nnode
      kss1(i,j) = (Fat_deg*Gc/Ls0 + 2.d0*KFreeEner)*shapef(i)*shapef(j)
      spf(i,j)  = ((Fat_deg*Gc*Ls0)*kss2(i,j) +
     *             kss1(i,j))*factor     
        enddo
      enddo


       RETURN
       END

!---------------------------------------------------------
!---------------------------------------------------------        
!---------------------------------------------------------    
      subroutine kStaticCondenS10(Kedd,Fed,KedEAS,KeEASEASinv,
     * KeEASd,FeEAS,ndim,nnode,nEAS)

        INTEGER i,j,nnode,ndim,nEAS
        
        REAL*8 Kedd(ndim*nnode,ndim*nnode),Fed(ndim*nnode),FeEAS(nEAS,1) 
        REAL*8 KedEAS(ndim*nnode,nEAS), KeEASEASinv(nEAS,nEAS), KeEASd(nEAS,ndim*nnode) 
        REAL*8 KeddAux(ndim*nnode,ndim*nnode),FedAux(ndim*nnode,1)
          
	  FedAux(:,1) = 0.d0
      KeddAux(:,:) = 0.d0    

      ! KeddAux =  KedEAS*KeEASEAS^-1*KeEASd
       call  kmatrixtriplS10(KeddAux,ndim*nnode,ndim*nnode,
     * KedEAS,ndim*nnode,nEAS,KeEASEASinv,nEAS,nEAS,
     * KeEASd,nEAS,ndim*nnode) 
c        This subroutine computing 
c        A = B*C*D
      ! FedAux  =  KedEAS*KeEASEAS^-1*FeEAS
       call  kmatrixtriplS10(FedAux,ndim*nnode,1,
     * KedEAS,ndim*nnode,nEAS,KeEASEASinv,nEAS,nEAS,
     * FeEAS,nEAS,1)  

		Kedd(:,:) = Kedd(:,:) - KeddAux(:,:) 
		Fed(:) = Fed(:) - FedAux(:,1)
      RETURN
      END

!---------------------------------------------------------
!---------------------------------------------------------         
!---------------------------------------------------------         
!---------------------------------------------------------         
!---------------------------------------------------------         
!---------------------------------------------------------         
!---------------------------------------------------------         
!---------------------------------------------------------
c*****************************************************************
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt,
     1 drplde,drpldt,stran,dstran,time,dtime,temp2,dtemp,predef,dpred,
     2 cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     3 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,kinc)

      use kvisual
        include 'aba_param.inc' !implicit real(a-h o-z)

      character*8 cmname
      dimension stress(ntens),statev(nstatv),ddsdde(ntens,ntens),
     1 ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens),
     2 time(2),predef(1),dpred(1),props(nprops),coords(3),drot(3,3),
     3 dfgrd0(3,3),dfgrd1(3,3),jstep(4)

		dimension knpt_array(npt)
		
      ddsdde=0.0d0
	  !noel=
	   knpt_array= (/ 5, 6, 7, 8, 1, 2, 3, 4/) !it is important to include comma
	  knpt=  knpt_array(npt)
	  
      noffset=noel-100000
      statev(1:nstatv-1)=USRVAR(noffset,1:nstatv-1,knpt)
      statev(nstatv)=Ac
      return
      end         
!---------------------------------------------------------         
!---------------------------------------------------------  