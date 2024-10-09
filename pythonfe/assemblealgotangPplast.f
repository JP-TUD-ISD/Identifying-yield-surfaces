! assembles material tangents to total tangent, would include total derivative of FP
      subroutine assemblealgotangPplast(dPdf, dPdfp, totFp, tang)

      implicit none

      integer i, j, k, l, m, n
      real (kind=8) dPdf(3,3,3,3),dPdfp(3,3,3,3), totFP(3,3,3,3)
      real (kind=8) tang(3,3,3,3)

      do i = 1,3
            do j = 1,3
                  do k = 1,3
                        do l = 1,3
                              tang(i,j,k,l) = dPdf(i,j,k,l)
                              do m = 1,3
                                    do n = 1,3
      tang(i,j,k,l) = tang(i,j,k,l) + dPdfp(i,j,m,n) * totFP(m,n,k,l)
                                    enddo ! n
                              enddo ! m
                        enddo ! l
                  enddo ! k
            enddo ! j
      enddo ! i

      end subroutine
