!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!   Single layer representation for Laplace equation in 2d
!
!  PDE:
!    \Delta u = 0
!
!  Representation:
!    u = S_{0}[\sigma]
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  User callable routines:
!    - lap2d_slp: Given a density sigma this routine returns the 
!        solution u
!
!    - lap2d_slp_eval: Given \sigma this routine evaluates the solution 
!        at a collection of targets
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  Advanced interfaces:
!    - getnearquad_lap2d_slp: compute the near quadrature correction
!        for constructing the on-surface Laplace single layer operator 
!        in 2d using user-provided near-information prescribed
!        in row-sparse compressed format
!
!    - getnearquad_lap2d_slp_eval: compute the near quadrature
!        correction for target points which can be either on-surface or
!        off-surface with user-provided near-information prescribed in
!        row-sparse compressed format
!

      subroutine getnearquad_lap2d_slp(npatches, norders, &
       ixyzs, iptype, npts, srccoefs, srcvals, eps, iquadtype, nnz, & 
       row_ptr, col_ind, iquad, rfac0, nquad, wnear)
!
!  This subroutine generates the near field quadrature
!  for the representation:
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
!        koornwinder expansion coefficients of xyz, dxyz/du,
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
!        pair. 
!
!  Output arguments
!    - wnear: real *8(nquad)
!        The desired near field quadrature
!               
!

      implicit none 
      integer, intent(in) :: npatches, norders(npatches)
      integer, intent(in) :: ixyzs(npatches+1), iptype(npatches)
      integer, intent(in) :: npts
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      real *8, intent(in) :: eps
      integer, intent(in) :: iquadtype, nnz
      integer, intent(in) :: row_ptr(npts+1), col_ind(nnz), iquad(nnz+1)
      real *8, intent(in) :: rfac0
      integer, intent(in) :: nquad
      
      real *8, intent(out) :: wnear(nquad)

      integer ndtarg, ntarg, ndd, ndz, ndi
      real *8 dpars(1)
      complex *16 zpars(1)
      integer ipars(1)

      real *8, allocatable :: uvs_targ(:,:)
      integer, allocatable :: ipatch_id(:)

      integer i
      integer ipv

      procedure (), pointer :: fker
      external l2d_slp


      ndz=0
      ndd=0
      ndi=0

      ndtarg = 12
      ntarg = npts

      allocate(ipatch_id(npts),uvs_targ(2,npts))
!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo
!$OMP END PARALLEL DO      

      call get_patch_id_uvs(npatches, norders, ixyzs, iptype, npts, &
        ipatch_id, uvs_targ)

      ipv=0
      fker => l2d_slp 
      call dgetnearquad_ggq_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, srcvals, &
        ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear)

      return
      end subroutine getnearquad_lap2d_slp
!
!
!
!
