! calculates inverse of a 3x3 matrix and determinant
      subroutine invert3new(a, deta)

      implicit   none

      integer    i, j
      real*8     a(3,3), ai(3,3), deta, deti

      ai(1,1) = a(2,2)*a(3,3) - a(3,2)*a(2,3)
      ai(1,2) = a(3,2)*a(1,3) - a(1,2)*a(3,3)
      ai(1,3) = a(1,2)*a(2,3) - a(2,2)*a(1,3)

      deta = a(1,1)*ai(1,1) + a(2,1)*ai(1,2) + a(3,1)*ai(1,3)

      deti    = 1.d0/deta
      ai(1,1) = ai(1,1)*deti
      ai(1,2) = ai(1,2)*deti
      ai(1,3) = ai(1,3)*deti

      ai(2,1) = (a(2,3)*a(3,1) - a(3,3)*a(2,1))*deti
      ai(2,2) = (a(3,3)*a(1,1) - a(1,3)*a(3,1))*deti
      ai(2,3) = (a(1,3)*a(2,1) - a(2,3)*a(1,1))*deti

      ai(3,1) = (a(2,1)*a(3,2) - a(3,1)*a(2,2))*deti
      ai(3,2) = (a(3,1)*a(1,2) - a(1,1)*a(3,2))*deti
      ai(3,3) = (a(1,1)*a(2,2) - a(2,1)*a(1,2))*deti

      do j = 1,3
        do i = 1,3
          a(i,j) = ai(i,j)
        enddo ! i
      enddo ! j

      end
