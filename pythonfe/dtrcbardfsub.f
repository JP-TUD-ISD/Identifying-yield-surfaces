! calculates derivative of trace of isochoric part of right Cauchy Green tensor with respect to deformation gradient
      subroutine dtrcbardfsub(j23, f, finv, trcbar, res)

      implicit none

      integer i,j

      real*8 j23, f(3,3), finv(3,3), trcbar, res(3,3)

      do i = 1,3
            do j = 1,3
                  res(i,j) = 2.d0 * (f(i,j) * j23 - trcbar/3.d0 * 
     &            finv(j,i))
            enddo ! j
      enddo ! i

      end subroutine
