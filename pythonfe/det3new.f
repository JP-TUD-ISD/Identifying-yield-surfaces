! calculates determinant of 3x3 matrix
      subroutine det3new(a, deta)

      implicit none

      real (kind=8) a(3,3), deta, ai(3,3)

      ai(1,1) = a(2,2)*a(3,3) - a(3,2)*a(2,3)
      ai(1,2) = a(3,2)*a(1,3) - a(1,2)*a(3,3)
      ai(1,3) = a(1,2)*a(2,3) - a(2,2)*a(1,3)

      deta = a(1,1)*ai(1,1) + a(2,1)*ai(1,2) + a(3,1)*ai(1,3)

      end subroutine
