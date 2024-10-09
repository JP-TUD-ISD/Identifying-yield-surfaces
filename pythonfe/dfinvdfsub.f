! calculates derivative of F^-1 with respect to F (deformation gradient)
      subroutine dfinvdfsub(finv, res)

      implicit none

      integer i,j,k,l

      real*8 finv(3,3), res(3,3,3,3)

      do i = 1,3
            do j = 1,3
                  do k = 1,3
                        do l = 1,3
      res(i,j,k,l) = - finv(i,k) * finv(l,j)
                        enddo ! l
                  enddo ! k
            enddo ! j
      enddo ! i

      end subroutine
