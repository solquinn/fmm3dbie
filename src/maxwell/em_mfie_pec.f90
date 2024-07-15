!
!  Magnetic Field integral equation (MFIE)
!  for scattering from perfect electric conductors,
!  with routines for computing the charge density
!
!  PDE:
!    \nabla \times E =  ik H
!    \nabla \cdot  E =     0
!    \nabla \times H = -ik E
!    \nabla \cdot  H =     0
!
!  Boundary conditions
!    n \times (E + E_in) = 0
!    n \cdot  (E + E_in) = \rho
!
!    n \times (H + H_in) = J
!    n \cdot  (H + H_in) = 0
!  
!
!  Representation:
!    H = \nabla \times S_{k} [J]
!    E = ik S_{k} [J] - \nabla S_{k} [\rho]
!
!  Integral equations solve for [J] and obtained by imposing
!
!  n \times (H + H_in) = J
!
!  and is given by
!  J/2 - n \times \nabla S_{k} [J] =  n \times H_in  ----  (1)
!
!  Once J is computed $\rho$ (if requested) is obtained by imposing
!
!  n \cdot  (E + E_in) = \rho
!  \rho/2 + S_{k}'[\rho] = n \cdot E_in + ik n \cdot S_{k}[J] 
!
!
!  Like all other maxwell formulations, the vector densties and equation
!  are expanded in a local orthonormal basis given by
!  \partial_{u} xyz, and \bn \times \partial_{u} xyz appropriately 
!  normalized. Let $\ru$ and $\rv$ denote those unit tangent vectors, 
!  then the currents are given by
!
!  so J = j_{u} \ru + j_{v} \rv  
!
!  and equation (1) is imposed by dotting it with \ru and \rv
!
!  User callable routines:
!    - em_mfie_pec_solver: Given incident electric and magnetic 
!        fields E_in, H_in, helmholtz wave number k,
!        and parameter \alpha, this routine returns the current J, and
!        charge \rho 
!
!    - em_mfie_pec_eval: Given J, and \rho,
!        evaluate the electric and magnetic fields at a collection of targets
!
!  Advanced interfaces:
!    - getnearquad_em_mfie_pec: Computes the quadrature
!        correction for constructing the on-surface integral equation
!        with user-provided near-information prescribed in 
!        row-sparse compressed format
!
!    - lpcomp_em_mfie_pec_addsub: Apply the principal value part
!        of the integral equation on surface. On input, user
!        provides precomputed near quadrature in row-sparse
!        compressed format, and oversampling information for 
!        handling the far part of the computation
!
!    - em_mfie_pec_solver_guru: Guru solver routine, 
!        where user is responsible for providing precomputed
!        near quadrature information in row-sparse compressed format
!        and oversampling surface information for the far-part
!
!    - getnearquad_em_mfie_pec_eval: Generate quadrature
!        for the post-processor, where targets can be in the volume
!        or on surface
!
!
!    - em_mfie_pec_eval_addsub: Compute the solution
!        E, H at a collection of targets, given J and
!        \rho. On input, user provides precomputed
!        near quadrature in row-sparse compressed format, 
!        and oversampling surface information for handling the
!        far part
!   



      subroutine getnearquad_em_mfie_pec(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, eps, zpars, &
        iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, &
        wnear)
!
!  This subroutine generates the near field quadrature
!  for the magnetic field integral equation.
!
!  The quadrature is computed by the following strategy
!  targets within a sphere of radius rfac0*rs
!  of a chunk centroid is handled using adaptive integration
!  where rs is the radius of the bounding sphere
!  for the patch
!  
!  All other targets in the near field are handled via
!  oversampled quadrature
!
!  The recommended parameter for rfac0 is 1.25d0
!
!  NOTES:
!  This subroutine returns 9 kernels
!  1) \ru \cdot (-n \times \nabla \times S_{k})[j_{u} \ru]
!  2) \ru \cdot (-n \times \nabla \times S_{k})[j_{v} \rv]
!  3) \rv \cdot (-n \times \nabla \times S_{k})[j_{u} \ru]
!  4) \rv \cdot (-n \times \nabla \times S_{k})[j_{v} \rv]
!   
!  Input arguments:
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization on each patch 
!    - ixyzs: integer(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!        iptype = 11, quadrangular patch discretized with GL nodes
!        iptype = 12, quadrangular patch discretized with Chebyshev 
!                     nodes
!    - npts: integer
!        total number of discretization points on the boundary
!    - srccoefs: real *8 (9,npts)
!        basis expansion coefficients of xyz, dxyz/du,
!        and dxyz/dv on each patch. 
!        For each point 
!          * srccoefs(1:3,i) is xyz info
!          * srccoefs(4:6,i) is dxyz/du info
!          * srccoefs(7:9,i) is dxyz/dv info
!    - srcvals: real *8 (12,npts)
!        xyz(u,v) and derivative info sampled at the 
!        discretization nodes on the surface
!          * srcvals(1:3,i) - xyz info
!          * srcvals(4:6,i) - dxyz/du info
!          * srcvals(7:9,i) - dxyz/dv info
!          * srcvals(10:12,i) - normals info
!    - eps: real *8
!        precision requested
!    - zpars: complex *16 (2)
!        kernel parameters
!          * zpars(1) = k 
!          * zpars(2) = alpha
!    - iquadtype: integer
!        quadrature type
!          * iquadtype = 1, use ggq for self + adaptive integration
!            for rest
!    - nnz: integer
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer(npts+1)
!        row_ptr(i) is the pointer
!        to col_ind array where list of relevant source patches
!        for target i start
!    - col_ind: integer (nnz)
!        list of source patches relevant for all targets, sorted
!        by the target number
!    - iquad: integer(nnz+1)
!        location in wnear_ij array where quadrature for col_ind(i)
!        starts for a single kernel. In this case the different kernels
!        are matrix entries are located at (m-1)*nquad+iquad(i), where
!        m is the kernel number
!    - rfac0: integer
!        radius parameter for near field
!    - nquad: integer
!        number of near field entries corresponding to each source target
!        pair. The size of wnear is 9*nquad since there are 9 kernels
!        per source target pair
!
!    output
!      wnear - complex *16 (4, nquad) 
!        Near field quadrature corrections for mfie
!        * wnear(1,:) stores corrections for
!          \ru \cdot (-n \times \nabla \times S_{k})[j_{u} \ru]
!        * wnear(2,:) stores corrections for
!          \ru \cdot (-n \times \nabla \times S_{k})[j_{v} \rv]
!        * wnear(3,:) stores corrections for
!          \rv \cdot (-n \times \nabla \times S_{k})[j_{u} \ru]
!        * wnear(4,:) stores corrections for
!          \rv \cdot (-n \times \nabla \times S_{k})[j_{v} \rv]
!

      implicit none 
      integer, intent(in) :: npatches, norders(npatches), npts, nquad
      integer, intent(in) :: ixyzs(npatches+1), iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts), eps
      real *8, intent(in) :: rfac0
      
      integer, intent(in) :: iquadtype, nnz
      integer, intent(in) :: row_ptr(npts+1), col_ind(nnz)
      integer, intent(in) :: iquad(nnz+1)
      
      complex *16, intent(in) :: zpars(2)
      complex *16, intent(out) :: wnear(4, nquad)
      
      integer ndtarg
      
      integer ipars(2)
      real *8 dpars(1)
      
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_src(:,:)
      complex *16, allocatable :: wneartmp(:)


      complex *16 alpha, beta
      integer i, j, ndi, ndd, ndz, iuse, idim, l

      integer ipv

      procedure (), pointer :: fker
      external  fker_em_mfie_pec

      ndz=1
      ndd=1
      ndi=2
      ipv=1
      ndtarg = 12

      fker =>  fker_em_mfie_pec

      allocate(ipatch_id(npts), uvs_src(2,npts))
      !$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        ipatch_id(i) = -1
        uvs_src(1,i) = 0
        uvs_src(2,i) = 0
      enddo
!$OMP END PARALLEL DO      

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, idim)
      do i=1,nquad
        do idim = 1,4
          wnear(idim,i) = 0
        enddo
      enddo
!$OMP END PARALLEL DO

      allocate(wneartmp(nquad))
      call get_patch_id_uvs(npatches, norders, ixyzs, iptype, npts, &
        ipatch_id, uvs_src)

      ipv = 1
      do j=1,2
        do i=1,2
          ipars(1) = j
          ipars(2) = i
          iuse = (j-1)*2 + i
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(l)
          do l=1,nquad
            wneartmp(l) = 0
          enddo
!$OMP END PARALLEL DO                   
          call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
            iptype, npts, srccoefs, srcvals, ndtarg, npts, srcvals, &
            ipatch_id, uvs_src, eps, ipv, fker, ndd, dpars, ndz, &
            zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, &
            nquad, wneartmp)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(l)
          do l=1,nquad
            wnear(iuse,l) = wneartmp(l)
          enddo
!$OMP END PARALLEL DO  
        enddo
      enddo
      
      return
      end subroutine getnearquad_em_mfie_pec
!
!
!
!
!

      subroutine lpcomp_em_mfie_pec_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, &
        novers, nptso, ixyzso, srcover, whtsover, lwork, work, &
        ndim, sigma, pot)
!
!
!  This subroutine evaluates the layer potential for
!  the boundary integral equation:
!
!  J/2 - n \times \nabla S_{k} [J] =  n \times H_in  
!
!  where the near field is precomputed and stored
!  in the row sparse compressed format.
!
!  Note the 4\pi scaling is NOT included as the FMM output was scaled
!  appropriately
!
!  Note: the identity J/2 is not included as the gmres takes care of
!  that
!
!  The fmm is used to accelerate the far-field and 
!  near-field interactions are handled via precomputed quadrature
!
!  Using add and subtract - no need to call tree and set fmm 
!  parameters can directly call existing fmm library
!
!  Input arguments:
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization on each patch 
!    - ixyzs: integer(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!        iptype = 11, quadrangular patch discretized with GL nodes
!        iptype = 12, quadrangular patch discretized with Chebyshev nodes
!    - npts: integer
!        total number of discretization points on the boundary
!    - srccoefs: real *8 (9,npts)
!        basis expansion coefficients of xyz, dxyz/du,
!        and dxyz/dv on each patch. 
!        For each point 
!          * srccoefs(1:3,i) is xyz info
!          * srccoefs(4:6,i) is dxyz/du info
!          * srccoefs(7:9,i) is dxyz/dv info
!    - srcvals: real *8 (12,npts)
!        xyz(u,v) and derivative info sampled at the 
!        discretization nodes on the surface
!          * srcvals(1:3,i) - xyz info
!          * srcvals(4:6,i) - dxyz/du info
!          * srcvals(7:9,i) - dxyz/dv info
!          * srcvals(10:12,i) - normals info
!    - eps: real *8
!        precision requested
!    - ndd: integer
!        number of real parameters defining the kernel/
!        integral representation (unused in this routine)
!    - dpars: real *8(ndd)
!        real parameters defining the kernel/
!        integral represnetation (unused in this routine)
!    - ndz: integer
!        number of complex parameters defining the kernel/
!        integral representation, must be >=1 
!    - zpars: complex *16(ndz)
!        complex parameters defining the kernel/
!        integral represnetation.
!        * zpars(1) = k 
!    - ndi: integer
!        number of integer parameters defining the kernel/
!        integral representation (unused in this routine)
!    - ipars: integer(ndi)
!        integer parameters defining the kernel/
!        integral represnetation (unused in this routine)
!    - nnz: integer
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer(npts+1)
!        row_ptr(i) is the pointer
!        to col_ind array where list of relevant source patches
!        for target i start
!    - col_ind: integer (nnz)
!        list of source patches relevant for all targets, sorted
!        by the target number
!    - iquad: integer(nnz+1)
!        location in wnear_ij array where quadrature for col_ind(i)
!        starts for a single kernel. In this case the different kernels
!        are matrix entries are located at (m-1)*nquad+iquad(i), where
!        m is the kernel number
!    - nquad: integer
!        number of near field entries corresponding to each source target
!        pair
!    - nker: integer
!        number of kernels in quadrature correction, must be 4
!    - wnear: complex *16(nker, nquad)
!        Near field quadrature corrections for mfie
!        * wnear(1,:) stores corrections for
!          \ru \cdot (-n \times \nabla \times S_{k})[j_{u} \ru]
!        * wnear(2,:) stores corrections for
!          \ru \cdot (-n \times \nabla \times S_{k})[j_{v} \rv]
!        * wnear(3,:) stores corrections for
!          \rv \cdot (-n \times \nabla \times S_{k})[j_{u} \ru]
!        * wnear(4,:) stores corrections for
!          \rv \cdot (-n \times \nabla \times S_{k})[j_{v} \rv]
!    - novers: integer(npatches)
!        order of discretization for oversampled sources and
!        density
!    - ixyzso: integer(npatches+1)
!        ixyzso(i) denotes the starting location in srcover,
!        corresponding to patch i
!    - nptso: integer
!        total number of oversampled points
!    - srcover: real *8 (12,nptso)
!        oversampled set of source information
!    - whtsover: real *8 (nptso)
!        smooth quadrature weights at oversampled nodes
!    - lwork: integer
!        size of work array (unused in this routine)
!    - work: real *8(lwork)
!        work array (unused in this routine)
!    - ndim: integer
!        number of densities (must be 2)
!    - sigma: complex *16(ndim, npts)
!        * sigma(1,:) is the ru component of J 
!        * sigma(2,:) is the rv component of J

!
!  Output arguments
!    - pot: complex *16 (ndim, npts)
!        The application of the integral representation 
!        (excluding the identity term)
!        * pot(1,:) is the ru component of n x H 
!        * pot(2,:) is the rv component of n x H 
!

      implicit none
      integer, intent(in) :: npatches, npts
      integer, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      real *8, intent(in) :: eps
      integer, intent(in) :: ndd, ndi, ndz
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zpars(ndz)
      integer, intent(in) :: ipars(ndi)
      integer, intent(in) :: nnz, nquad
      integer, intent(in) :: row_ptr(npts+1), col_ind(nnz)
      integer, intent(in) :: iquad(nnz+1)
      
      integer, intent(in) :: ndim
      complex *16, intent(in) :: sigma(ndim,npts)
      
      integer, intent(in) :: nptso
      integer, intent(in) :: novers(npatches+1)
      integer, intent(in) :: ixyzso(npatches+1)
      real *8, intent(in) :: srcover(12,nptso), whtsover(nptso)
      
      integer, intent(in) :: nker
      complex *16 wnear(nker, nquad)

      integer, intent(in) :: lwork
      real *8, intent(in) :: work(lwork)

      complex *16, intent(out) :: pot(ndim,npts)

!
!  temporary variables
!
      integer nduse 
      real *8, allocatable :: ru(:,:), rv(:,:)
      real *8, allocatable :: ruover(:,:), rvover(:,:)

      complex *16, allocatable :: zjvec(:,:), charges(:,:)
      complex *16, allocatable :: zpot(:,:), zgrad(:,:,:)
      complex *16, allocatable :: sigmaover(:,:)
      integer nmax
      real *8 thresh,ra
      real *8, allocatable :: sources(:,:)

      complex *16 alpha, zk 

      real *8, allocatable :: srctmp(:,:), srctmp2(:,:)
      complex *16, allocatable :: ctmp(:,:)
      complex *16 zpottmp(3), zgradtmp(3,3), pottmp(2)

      integer ns, nt, npols
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      integer i,j,jpatch,jquadstart,jstart

      real *8 timeinfo(10),t1,t2
!$    real *8 omp_get_wtime      

      real *8 over4pi
      integer nss,ii,l,npover, m, n, iind, ier

      integer nd,ntarg0

      real *8 ttot,done,pi

      parameter (nd=1,ntarg0=1)
      data over4pi/0.07957747154594767d0/


      ns = nptso
      done = 1
      pi = atan(done)*4
           
      ifpgh = 0
      ifpghtarg = 1
      nduse = 3
      allocate(sources(3,ns), srctmp(3,npts))
      allocate(charges(nduse,ns), zjvec(3,ns))
      allocate(zpot(nduse,npts), zgrad(nduse,3,npts))
      allocate(sigmaover(2,ns))
      allocate(ru(3,npts), rv(3,npts))
      allocate(ruover(3,nptso), rvover(3,nptso))

      call orthonormalize_all(srcvals(4:6,:), srcvals(10:12,:), ru, &
         rv, npts)
      call orthonormalize_all(srcover(4:6,:), srcover(10:12,:), &
         ruover, rvover, ns)

      call cpu_time(t1)
!$      t1=omp_get_wtime()      
      call get_near_corr_max(npts, row_ptr, nnz, col_ind, npatches, &
        ixyzso, nmax)
      call cpu_time(t2)
!$      t2=omp_get_wtime()     

      call cpu_time(t1)
!$      t1=omp_get_wtime()      

! 
!       oversample density
!

      call oversample_fun_surf(4, npatches, norders, ixyzs, iptype, & 
         npts, sigma, novers, ixyzso, ns, sigmaover)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,ns
         sources(1:3,i) = srcover(1:3,i)
         zjvec(1:3,i) = sigmaover(1,i)*ruover(1:3,i) + &
              sigmaover(2,i)*rvover(1:3,i)
         charges(1:3,i) = zjvec(1:3,i)*whtsover(i)*over4pi
      enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,npts
         srctmp(1:3,i) = srcvals(1:3,i)
      enddo
!$OMP END PARALLEL DO


!
!        compute threshold for ignoring local computation
!
      call get_fmm_thresh(12, ns, srcover, 12, npts, srcvals, thresh)
      call cpu_time(t2)
!$      t2=omp_get_wtime()      
      
      ier = 0
      zk = zpars(1)
      call hfmm3d_t_c_g_vec(nduse, eps, zk, ns, sources, &
        charges, npts, srctmp, zpot, zgrad, ier)
!
! Convert cartesian fmm info into components of integral equation
!

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) 
      do i=1,npts
         call get_mfie_inteq_comps_from_potgrad(zk, &
            srcvals(1,i), ru(1,i), rv(1,i), zpot(1,i), zgrad(1,1,i), &
            pot(1,i))
      enddo
!$OMP END PARALLEL DO

      
!
!       Add near field precomputed contribution
!
      call cpu_time(t1)
!$      t1 = omp_get_wtime()

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, jpatch, jquadstart) &
!$OMP PRIVATE(jstart, npols, l, m, n, iind)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            do m=1,2
              do n=1,2
                iind = (m-1)*2 + n
                pot(m,i) = pot(m,i) + wnear(iind,jquadstart+l-1)*  &
                             sigma(n,jstart+l-1)
              enddo
            enddo
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

      call cpu_time(t2)
!$      t2 = omp_get_wtime()
      timeinfo(2) = t2-t1

!
!     Remove near contribution of the FMM
!
      allocate(srctmp2(3,nmax), ctmp(nduse,nmax))
      call cpu_time(t1)
!$      t1 = omp_get_wtime()
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, nss, j, jpatch, l)  &
!$OMP PRIVATE(srctmp2, ctmp, zpottmp, zgradtmp, pottmp)
      do i=1,npts
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l=ixyzso(jpatch), ixyzso(jpatch+1)-1
            nss = nss + 1
            srctmp2(1:3,nss) = sources(1:3,l)            
            ctmp(1:nduse,nss) = charges(1:nduse,l)
          enddo
        enddo
        zpottmp(1:3) = 0
        zgradtmp(1:3,1:3) = 0
        pottmp(1:2) = 0
        call h3ddirectcg(nduse, zk, srctmp2, ctmp, nss, &
          srctmp(1,i), ntarg0, zpottmp, zgradtmp, thresh)

        call get_mfie_inteq_comps_from_potgrad(zk, &
            srcvals(1,i), ru(1,i), rv(1,i), zpottmp, zgradtmp, &
            pottmp)

        pot(1:2,i) = pot(1:2,i) - pottmp(1:2)
      enddo
!$OMP END PARALLEL DO      
      call cpu_time(t2)
!$      t2 = omp_get_wtime()     

      timeinfo(3) = t2-t1
      
      
      ttot = timeinfo(1) + timeinfo(2) +timeinfo(3)


      return
      end subroutine lpcomp_em_mfie_pec_addsub
!
!
!
!
!

      subroutine get_mfie_inteq_comps_from_potgrad(zk, srcvals, &
            ru, rv, zpot, zgrad, pot)
!
!  Given the potential and gradient of the vector helmholtz FMM which computes
!  the cartesian coordinates of S_{k}[J], this subroutine
!  returns the components of the mfie integral equation given by
!
!    -M_{k}[J]  -  ru component 
!    -M_{k}[J]  -  rv component 
!  
!  where M_{k} = n \times \nabla \times S_{k} is the magnetic field 
!  integral equation operator
!  
!  Input Arguments:
!    - zk: complex *16
!         Wavenumber (k) above
!    - srcvals: real *8(12)
!         target info, xyz coordinates, d/du (xyz), d/dv(xyz), and 
!         normals
!    - ru: real *8(3)
!         orthonormalized tangent vectors, d/du(xyz)/|d/du(xyz)|
!    - rv: real *8(3)
!         second orthonormalized tangent vector
!         ru \times normal
!    - zpot: complex *16(3)
!         the 4 components of the potential, 
!         S_{k}[J_{1}], S_{k}[J_{2}], S_{k}[J_{3}], S_{k}[\rho]
!    - zgrad: complex *16(3,3)
!         gradients of the potential above
!
!
!  Output arguments:
!    - pot: complex *16(2)
!        the components of the integral equation
!
!
!
      implicit none
      complex *16, intent(in) :: zk
      real *8, intent(in) :: srcvals(12), ru(3), rv(3)
      complex *16, intent(in) :: zpot(3), zgrad(3,3)
      complex *16, intent(out) :: pot(2)

      complex *16 zcurl(3), zvec(3)


      zcurl(1) = zgrad(3,2) - zgrad(2,3)
      zcurl(2) = zgrad(1,3) - zgrad(3,1)
      zcurl(3) = zgrad(2,1) - zgrad(1,2)

      call dzcross_prod3d(srcvals(10), zcurl, zvec)

      pot(1) = -(zvec(1)*ru(1) + zvec(2)*ru(2) + zvec(3)*ru(3))
      pot(2) = -(zvec(1)*rv(1) + zvec(2)*rv(2) + zvec(3)*rv(3))
         

      return
      end
!
!
!
!
!
!
      subroutine em_mfie_pec_solver(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpars, numit, einc, &
        hinc, eps_gmres, niter, errs, rres, zjvec, ifaugsolve, rho)
!
!
!  This subroutine solves the Scattering Maxwell p.e.c. problem.
!  The the equations are:
!
!    curlE = ikH; curl H = -ikE
!
!  Representation:
!
!    H = curlS_{k}[J]
!
!    E = ikS_{k}[J] - gradS_{k}[rho]
!
!  Boundary conditions:
!
!    (1)                       b.c.1
!  
! ifaugsolve = 1, then
!    (3)                       b.c.2
!
!    where:
!
!    n x H + n x H_inc = J     (1)
!
!    n x E + n x E_inc = 0     (2)
!
!    n · E + n · E_inc = rho   (3)
!
!    (div E)/ik = 0                 (4)
!
!  The incoming fields must be 'compatible' (solutions to Maxwell's 
!  equations in free space)
!
!  Boundary integral equation:
!
!    J/2 - M_{k}[J] = nxH_inc 
!    rho/2 + S_{k}'[rho] = n·E_inc + ik n \cdot S_{k}[J]
!     
!
!  The linear system is solved iteratively using GMRES
!  until a relative residual of eps_gmres is reached
!
!
!  input:
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization on each patch 
!    - ixyzs: integer(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!        iptype = 11, quadrangular patch discretized with GL nodes
!        iptype = 12, quadrangular patch discretized with Chebyshev nodes
!    - npts: integer
!        total number of discretization points on the boundary
!    - srccoefs: real *8 (9,npts)
!        basis expansion coefficients of xyz, dxyz/du,
!        and dxyz/dv on each patch. 
!        For each point 
!          * srccoefs(1:3,i) is xyz info
!          * srccoefs(4:6,i) is dxyz/du info
!          * srccoefs(7:9,i) is dxyz/dv info
!    - srcvals: real *8 (12,npts)
!        xyz(u,v) and derivative info sampled at the 
!        discretization nodes on the surface
!          * srcvals(1:3,i) - xyz info
!          * srcvals(4:6,i) - dxyz/du info
!          * srcvals(7:9,i) - dxyz/dv info
!          * srcvals(10:12,i) - normals info
!    - eps: real *8
!      precision requested for computing quadrature and fmm
!      tolerance
!    - zpars: complex *16 (1)
!      kernel parameters (Referring to formula (1))
!      zpars(1) = k 
!    - numit: integer
!      max number of gmres iterations        
!    - einc: complex *16(3, npts)
!      incident electric field
!    - hinc: complex *16(3, npts)
!      incident magnetic field
!    - eps_gmres: real *8
!      gmres tolerance requested
!    - ifaugsolve: integer
!      flag for doing an augmented solve for computing the charge
!      density. 
!      ifaugsolve = 1, charge density computed, else
!          charge density not computed
!
!    output
!    - niter: integer
!        number of gmres iterations required for relative residual
!        to converge to eps_gmres
!    - errs:  real *8 (numit+1)
!        relative residual as a function of iteration
!        number (only errs(1:niter) is populated))
!    - rres: real *8
!        relative residual for computed solution
!    - zjvec: complex *16(3,npts)
!        the induced surface current
!    - rho: complex *16(npts)
!        the surface charge density, if requested          
!

      implicit none
      integer, intent(in) :: npatches, npts
      integer, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      real *8, intent(in) :: eps, eps_gmres
      complex *16, intent(in) :: zpars(2)
      complex *16, intent(in) :: hinc(3,npts), einc(3,npts)
      integer, intent(in) :: numit
      integer, intent(in) :: ifaugsolve

      real *8, intent(out) :: errs(numit+1), rres
      integer, intent(out) :: niter
      complex *16, intent(out) :: zjvec(3,npts), rho(npts)

      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_src(:,:)
    

      integer nover, npolso, nptso
      integer nnz, nquad
      integer, allocatable :: row_ptr(:), col_ind(:), iquad(:)

      complex *16, allocatable :: wnear(:,:)

      real *8, allocatable :: srcover(:,:), wover(:)
      integer, allocatable :: ixyzso(:), novers(:)

      real *8, allocatable :: cms(:,:), rads(:), rad_near(:) 

      integer i, j, jpatch, jquadstart, jstart

      integer ndtarg, nker
      real *8 ttot
      real *8 rfac, rfac0
      integer iptype_avg, norder_avg
      integer ikerorder, iquadtype, npts_over

!
!
!       gmres variables
!
      complex *16 zid
      integer iind, k, l
!
!
! 
      allocate(uvs_src(2, npts), ipatch_id(npts))
!C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,npts
        ipatch_id(i) = -1
        uvs_src(1,i) = 0
        uvs_src(2,i) = 0
      enddo
!C$OMP END PARALLEL DO   

!
!    initialize patch_id and uv_targ for on surface targets
!
      call get_patch_id_uvs(npatches, norders, ixyzs, iptype, npts, &
        ipatch_id, uvs_src)

!
!
!        this might need fixing
!
      iptype_avg = floor(sum(iptype)/(npatches+0.0d0))
      norder_avg = floor(sum(norders)/(npatches+0.0d0))

      call get_rfacs(norder_avg, iptype_avg, rfac, rfac0)


      allocate(cms(3,npatches), rads(npatches), rad_near(npatches))

      call get_centroid_rads(npatches, norders, ixyzs, iptype, npts, & 
        srccoefs, cms, rads)

!$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,npatches
        rad_near(i) = rads(i)*rfac
      enddo
!$OMP END PARALLEL DO      

!
!    find near quadrature correction interactions
!
      print *, "entering find near mem"
      ndtarg = 12
      call findnearmem(cms, npatches, rad_near, ndtarg, srcvals, npts, &
        nnz)

      allocate(row_ptr(npts+1), col_ind(nnz))
      
      call findnear(cms, npatches, rad_near, ndtarg, srcvals, npts, &
        row_ptr, col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches, ixyzs, npts, nnz, row_ptr, col_ind, &
        iquad)

      ikerorder = 0

!
!    estimate oversampling for far-field, and oversample geometry
!

      allocate(novers(npatches), ixyzso(npatches+1))

      print *, "beginning far order estimation"

      call get_far_order(eps, npatches, norders, ixyzs, iptype, cms, &
        rads, npts, srccoefs, ndtarg, npts, srcvals, ikerorder, &
        zpars(1), nnz, row_ptr, col_ind, rfac, novers, ixyzso)

      npts_over = ixyzso(npatches+1)-1
      

      allocate(srcover(12,npts_over), wover(npts_over))

      call oversample_geom(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, srcvals, novers, ixyzso, npts_over, srcover)

      call get_qwts(npatches, novers, ixyzso, iptype, npts_over, &
        srcover, wover)


!
!   compute near quadrature correction
!
      nquad = iquad(nnz+1)-1
      print *, "nquad=",nquad
      allocate(wnear(4,nquad))
      
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j)      
      do i=1,nquad
        do j=1,4
          wnear(j,i)=0
        enddo
      enddo
!$OMP END PARALLEL DO          


      iquadtype = 1

      call getnearquad_em_mfie_pec(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, eps, zpars, &
        iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, &
        wnear)
!
!
!  CONTINUE FROM HERE:
!  First add near quadrature for Sk and Sk' if ifaugsolve = 1
!  and then update the guru routine
!
      
      print *, "done generating near quadrature, now starting gmres"
      nker = 4

      call em_mfie_pec_solver_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpars, numit, &
        einc, hinc, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, &
        novers, npts_over, ixyzso, srcover, wover, eps_gmres, niter, &
        errs, rres, zjvec, ifaugsolve, rho)

!
      return
      end subroutine em_mfie_pec_solver
!
!
!
!
!                        
      subroutine em_mfie_pec_solver_guru(npatches, norders, ixyzs, &
      iptype, npts, srccoefs, srcvals, eps, zpars, numit, &
      einc, hinc, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, &
      novers, nptso, ixyzso, srcover, whtsover, eps_gmres, niter, &
      errs, rres, zjvec, ifaugsolve, rho)
!
!
!  This subroutine solves the Scattering Maxwell p.e.c. problem.
!  The the equations are:
!
!    curlE = ikH; curlH = -ikE
!
!  Representation:
!
!    H = curlS_{k}[J]
!
!    E = ikS_{k}[J]-gradS_{k}[rho]
!
!  Boundary conditions:
!
!    (1)                       b.c.1
!
!    (3)                       b.c.2
!
!    where:
!
!    n x H + n x H_inc = J     (1)
!
!    n x E + n x E_inc = 0     (2)
!
!    n · E + n · E_inc = rho   (3)
!
!    (div E)/ik = 0                 (4)
!
!  The incoming fields must be 'compatible' (solutions to Maxwell's 
!  equations in free space)
!
!  Boundary integral equations:
!
!    J/2 - M_{k}[J] = nxH_inc 
!    rho/2 + S_{k}'[rho] = n·E_inc + ik n \cdot S_{k}[J]
!
!  The linear system is solved iteratively using GMRES
!  until a relative residual of eps_gmres is reached.
!
!  This is the guru solver routine which on input
!  takes in precomputed near-quadrature corrections, 
!  and oversampled surface information for the far part                            
!  
!
!  input:
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization on each patch 
!    - ixyzs: integer(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!        iptype = 11, quadrangular patch discretized with GL nodes
!        iptype = 12, quadrangular patch discretized with Chebyshev nodes
!    - npts: integer
!        total number of discretization points on the boundary
!    - srccoefs: real *8 (9,npts)
!        basis expansion coefficients of xyz, dxyz/du,
!        and dxyz/dv on each patch. 
!        For each point 
!          * srccoefs(1:3,i) is xyz info
!          * srccoefs(4:6,i) is dxyz/du info
!          * srccoefs(7:9,i) is dxyz/dv info
!    - srcvals: real *8 (12,npts)
!        xyz(u,v) and derivative info sampled at the 
!        discretization nodes on the surface
!          * srcvals(1:3,i) - xyz info
!          * srcvals(4:6,i) - dxyz/du info
!          * srcvals(7:9,i) - dxyz/dv info
!          * srcvals(10:12,i) - normals info
!    - eps: real *8
!      precision requested for computing quadrature and fmm
!      tolerance
!    - zpars: complex *16 (1)
!      kernel parameters (Referring to formula (1))
!      zpars(1) = k 
!    - numit: integer
!      max number of gmres iterations
!    - einc: complex *16(3, npts)
!      incident electric field              
!    - hinc: complex *16(3, npts)
!      incident magnetic field
!    - nnz: integer
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer(npts+1)
!        row_ptr(i) is the pointer
!        to col_ind array where list of relevant source patches
!        for target i start
!    - col_ind: integer (nnz)
!        list of source patches relevant for all targets, sorted
!        by the target number
!    - iquad: integer(nnz+1)
!        location in wnear_ij array where quadrature for col_ind(i)
!        starts for a single kernel. In this case the different kernels
!        are matrix entries are located at (m-1)*nquad+iquad(i), where
!        m is the kernel number
!    - nquad: integer
!        number of near field entries corresponding to each source target
!        pair
!    - nker: integer
!        number of kernels in quadrature correction
!    - wnear: real *8(nker, nquad)
!        precomputed quadrature corrections 
!        Near field quadrature corrections for mfie
!        * wnear(1,:) stores corrections for
!          \ru \cdot (-n \times \nabla \times S_{k})[j_{u} \ru]
!        * wnear(2,:) stores corrections for
!          \ru \cdot (-n \times \nabla \times S_{k})[j_{v} \rv]
!        * wnear(3,:) stores corrections for
!          \rv \cdot (-n \times \nabla \times S_{k})[j_{u} \ru]
!        * wnear(4,:) stores corrections for
!          \rv \cdot (-n \times \nabla \times S_{k})[j_{v} \rv]
!    - novers: integer(npatches)
!        order of discretization for oversampled sources and
!        density
!    - ixyzso: integer(npatches+1)
!        ixyzso(i) denotes the starting location in srcover,
!        corresponding to patch i
!    - nptso: integer
!        total number of oversampled points
!    - srcover: real *8 (12,nptso)
!        oversampled set of source information
!    - whtsover: real *8 (nptso)
!        smooth quadrature weights at oversampled nodes
!    - eps_gmres: real *8
!      gmres tolerance requested
!    - ifaugsolve: integer
!      flag for doing an augmented solve for computing the charge
!      density. 
!      ifaugsolve = 1, charge density computed, else
!          charge density not computed
!
!    output
!    - niter: integer
!        number of gmres iterations required for relative residual
!        to converge to eps_gmres
!    - errs:  real *8 (numit+1)
!        relative residual as a function of iteration
!        number (only errs(1:niter) is populated))
!    - rres: real *8
!        relative residual for computed solution
!    - zjvec: complex *16(3,npts)
!        the induced surface current
!    - rho: complex *16(npts)
!        the surface charge density          
!
      implicit none
      integer, intent(in) :: npatches, npts
      integer, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      real *8, intent(in) :: eps
      complex *16, intent(in) :: zpars(2)
      integer, intent(in) :: numit
      complex *16, intent(in) :: hinc(3,npts), einc(3,npts)
      integer, intent(in) :: nnz, nquad, nker
      integer, intent(in) :: row_ptr(npts+1), col_ind(nnz), iquad(nnz+1)
      complex *16, intent(in) :: wnear(nker,nquad)
      integer, intent(in) :: nptso
      integer, intent(in) :: novers(npatches), ixyzso(npatches+1)
      real *8, intent(in) :: srcover(12,nptso), whtsover(nptso)
      real *8, intent(in) :: eps_gmres
      integer, intent(in) :: ifaugsolve
      integer, intent(out) :: niter
      real *8, intent(out) :: errs(numit+1), rres
      complex *16, intent(out) :: zjvec(3,npts), rho(npts)
      
      integer ndd, ndi, ndz
      integer ipars
      real *8 dpars, work
      complex *16, allocatable :: rhs(:,:), soln(:,:)
      real *8, allocatable :: ru(:,:), rv(:,:), wts(:)
      complex *16 zvec(3), zvec2(3), zvec3(3), ztmp, zid

      complex *16, allocatable :: wnear_skj(:), wnear_aug(:)
      complex *16 zpars_em(7)

      integer i, j, k, l, ndim, idim, lwork, ndtarg, nker_em
      integer idensflag, ipotflag, ndim_s, ndim_p
      

      complex *16, allocatable :: sigma_em(:,:), skj(:,:)

      procedure (), pointer :: fker
      external lpcomp_em_mfie_pec_addsub

      allocate(rhs(2,npts), soln(2,npts))
      allocate(ru(3,npts), rv(3,npts), wts(npts))

      call orthonormalize_all(srcvals(4:6,:), srcvals(10:12,:), ru, &
         rv, npts)

!
!  Initialize rhs from e and h fields
!                 
   
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(zvec, zvec2, ztmp)
      do i=1,npts

!  Compute n \times hinc        
        call dzcross_prod3d(srcvals(10,i), hinc(1,i), zvec)
       
        rhs(1,i) = ru(1,i)*zvec(1) + ru(2,i)*zvec(2) + ru(3,i)*zvec(3)
        rhs(2,i) = rv(1,i)*zvec(1) + rv(2,i)*zvec(2) + rv(3,i)*zvec(3)
      
        soln(1:2,i) = 0
      enddo
!$OMP END PARALLEL DO    
      
      ndd = 0
      ndz = 1
      ndi = 0
      lwork = 0

      fker => lpcomp_em_mfie_pec_addsub

      ndim = 2

      call get_qwts(npatches, norders, ixyzs, iptype, npts, &
      srcvals, wts)
!
      zid = 0.5d0
      call zgmres_guru(npatches, norders, ixyzs, &
            iptype, npts, srccoefs, srcvals, wts, &
            eps, ndd, dpars, ndz, zpars, ndi, ipars, &
            nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
            nptso, ixyzso, srcover, whtsover, lwork, work, &
            ndim, fker, zid, rhs, numit, eps_gmres, niter, errs, &
            rres, soln)
!
!  Now construct current and charge densities from solution
!                         

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,npts
        zjvec(1:3,i) = soln(1,i)*ru(1:3,i) + soln(2,i)*rv(1:3,i)
        rho(i) = 0
      enddo
!$OMP END PARALLEL DO            

!
!  Now solve for the density \rho if requested
!
      if (ifaugsolve.eq.1) then
!
!  First compute S_{k}[J] which is required data for the augmented
!  solve. 
!
!
        allocate(wnear_skj(nquad))
        call getnearquad_helm_comb_dir(npatches, norders, ixyzs, iptype, &
          npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, &
          row_ptr, col_ind, iquad, rfac0, nquad, wnear_skj)
        
        zpars_em(1) = zpars(1)
        zpars_em(2) = 1.0d0
        zpars_em(3) = 0.0d0
        zpars_em(4) = 0.0d0
        zpars_em(5) = 0.0d0
        zpars_em(6) = 0.0d0
        zpars_em(7) = 0.0d0

        ndz = 7

        ndtarg = 12

        nker_em = 1
        idensflag = 1
        ipotflag = 1

        allocate(sigma_em(8,npts), skj(3,npts))

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)        
        do i = 1,npts
          sigma_em(1,i) = zjvec(1,i)
          sigma_em(2,i) = zjvec(2,i)
          sigma_em(3,i) = zjvec(3,i)
          sigma_em(4,i) = 0
          sigma_em(5,i) = 0
          sigma_em(6,i) = 0
          sigma_em(7,i) = 0
          sigma_em(8,i) = 0

          skj(1,i) = 0
          skj(2,i) = 0
          skj(3,i) = 0
        enddo
!$OMP END PARALLEL DO

        ndim_s = 8
        ndim_p = 3
        call em_debye_eval_addsub(npatches, norders, ixyzs, iptype, &
          npts, srccoefs, srcvals, ndtarg, npts, srcvals, eps, ndd, &
          dpars, ndz, zpars_em, ndi, ipars, nnz, row_ptr, col_ind, &
          iquad, nquad, nker_em, wnear_skj, novers, nptso, ixyzso, &
          srcover, whtsover, lwork, work, idensflag, sigma_em, &
          ipotflag, ndim_p, skj)

      endif

      return
      end
!
!
!
!
!


      subroutine fker_em_mfie_pec(srcinfo, ndt, targinfo, ndd, dpars, &
       ndz, zpars, ndi, ipars, E_val)
!
! this function provides the near field kernel that will use 
! zgetnearquad_ggq_guru through getnearquad_em_mfie_pec
!
      implicit none
      integer, intent(in) :: ndt, ndd, ndz, ndi
      real *8, intent(in) :: srcinfo(12)
      real *8, intent(in) :: targinfo(ndt)
   	  integer, intent(in) :: ipars(ndi)
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zpars(ndz)
      complex *16, intent(out) :: E_val

! List of local variables

      complex *16 E_mat(3,3)
      real *8 ru_s(3), rv_s(3), n_s(3)
      real *8 ru_t(3), rv_t(3), n_t(3)
      real *8 du(3), dv(3),sour(3), targ(3)
      real *8 r, dr(3), aux_real

      complex *16 nxcurlSka(2,2), nxnxSkb(2,2)
      complex *16 ngradSklambda, nSkb(1,2)

      real *8 xprod_aux1(3), xprod_aux2(3), xprod_aux3(3)
      real *8 xprod_aux4(3)
      complex *16 R1, ima,my_exp, zk, alpha
      real *8 over4pi
      
      data over4pi/0.07957747154594767d0/
      data ima/(0.0d0, 1.0d0)/

      zk=zpars(1)
      alpha=zpars(2)
    
      sour(1)=srcinfo(1)
      sour(2)=srcinfo(2)
      sour(3)=srcinfo(3)

      n_s(1)=srcinfo(10)
      n_s(2)=srcinfo(11)
      n_s(3)=srcinfo(12)

      targ(1)=targinfo(1)
      targ(2)=targinfo(2)
      targ(3)=targinfo(3)

      n_t(1)=targinfo(10)
      n_t(2)=targinfo(11)
      n_t(3)=targinfo(12)

      dr(1)=targ(1)-sour(1)
      dr(2)=targ(2)-sour(2)
      dr(3)=targ(3)-sour(3)
      r=sqrt((dr(1))**2+(dr(2))**2+(dr(3))**2)
      
      R1=(ima*zk*r-1.0d0)/r**3*exp(ima*zk*r)*over4pi
      my_exp=exp(ima*zk*r)*over4pi

      call orthonormalize(srcinfo(4:6),n_s,ru_s,rv_s)
      call orthonormalize(targinfo(4:6),n_t,ru_t,rv_t)

      call get_nxcurlSka(ru_s, rv_s, n_s, ru_t, rv_t, n_t, dr, R1, &
        my_exp, r, nxcurlSka)

      E_mat(1,1)=-nxcurlSka(1,1)
      E_mat(1,2)=-nxcurlSka(1,2)
      E_mat(2,1)=-nxcurlSka(2,1)
      E_mat(2,2)=-nxcurlSka(2,2)

      call get_nSkb(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,my_exp,r,nSkb)
      E_mat(3,1) = ima*zk*nSkb(1,1)
      E_mat(3,2) = ima*zk*nSkb(1,2)

      E_val=E_mat(ipars(1),ipars(2))

      return
      end subroutine fker_em_mfie_pec




