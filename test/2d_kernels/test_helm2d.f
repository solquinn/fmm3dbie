      implicit real *8 (a-h,o-z)
      real *8 src(3), targ(12)
      complex *16 ima,z0,zk
      complex *16 val1,val2
      data ima /(0,1)/

      zk = 3.3d0

      src(1) = 0
      src(2) = 0 
      src(3) = 0

      targ(1) = 0.5d0
      targ(2) = -0.25d0
      targ(3) = 0
      targ(10) = 0.3d0
      targ(11) = 1.1d0

      call h2d_slp(src,12,targ,0,0,1,zk,0,0,val1)
      call h2d_sprime(src,12,targ,0,0,1,zk,0,0,val2)
      
      open(unit=33,file='print_test_helm2d.txt')
      write(33,'(2x,e22.16)'),val1
      write(33,'(2x,e22.16)'),val2
      close(33)

      stop
      end
