      subroutine sweight(lint,s)

      implicit  none

      integer   i,j,k,lint, ig(4),jg(4)
      real*8    g, s(4,*), sw(2,5),sqt13

      data      ig/-1,1,1,-1/,jg/-1,-1,1,1/

      sqt13 = dsqrt(1.d0/3.d0)

        lint = 8
        g    = sqt13
        do i = 1,4
          s(1,i)   = ig(i)*g
          s(1,i+4) = s(1,i)
          s(2,i)   = jg(i)*g
          s(2,i+4) = s(2,i)
          s(3,i)   =  g
          s(3,i+4) = -g
          s(4,i)   = 1.d0
          s(4,i+4) = 1.d0
        end do ! i

      end
