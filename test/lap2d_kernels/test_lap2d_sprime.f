      implicit real *8 (a-h,o-z)
      real *8 src(3), targ(12)
      complex *16 ima,z0,refval
      real *8 val
      data ima /(0,1)/

      src(1) = 0
      src(2) = 0
      src(3) = 0

      targ(1) = 0.5
      targ(2) = -0.25
      targ(3) = 0

      targ(10) = 1
      targ(11) = 0

      call l2d_sprime(src,12,targ,0,0,0,0,0,0,val)

      open(unit=33,file='print_test_lap2d.txt')
      write(33,'(2x,e22.16)'),val
      close(33)

      stop
      end
