! calculates first Piola Kirchhoff stresses and Mandel stresses in the intermediate configuration
! needs derivative of strain energy with respect to F^e
      subroutine calcPandMandel(dpsidfe,dpsidfefe, fe, fp, fpinv, Piola, 
     & mandel, dPioladfp, dfedfp, dfedf, dfedfdfp, dMdfp, dPdf, dMdf)

      implicit none

      integer i,j,k,l, m, n, o, p
      real (kind=8) dpsidfe(3,3), dpsidfefe(3,3,3,3), fe(3,3)
      real (kind=8) fp(3,3), Piola(3,3), Mandel(3,3)
      real (kind=8) dPioladfp(3,3,3,3), dfedfp(3,3,3,3), id(3,3)
      real (kind=8) dfedf(3,3,3,3), fpinv(3,3), dfedfdfp(3,3,3,3,3,3)
      real (kind=8) dMdfp(3,3,3,3), dPdf(3,3,3,3), dMdf(3,3,3,3)

      data id/ 1.d0 , 0.d0 , 0.d0 ,
     *         0.d0 , 1.d0 , 0.d0 ,
     *         0.d0 , 0.d0 , 1.d0 /

      do i = 1,3
            do j = 1,3
                  do k = 1,3
                        do l = 1,3
      dfedfp(i,j,k,l) = - fe(i,k) * fpinv(l,j)

      dfedf(i,j,k,l) = id(i,k) * fpinv(l,j)
                              do m = 1,3
                                    do n = 1,3
      dfedfdfp(i,j,k,l,m,n) = - id(i,k) * fpinv(l,m) * fpinv(n,j)
                                    enddo ! n
                              enddo ! m
                        enddo ! l
                  enddo ! k
            enddo ! j
      enddo ! i

      do i = 1,3
            do j = 1,3
                  Piola(i,j) = 0.d0
                  do k = 1,3
                        do l = 1,3
      Piola(i,j) = Piola(i,j) + dpsidfe(k,l) * dfedf(k,l,i,j)
                        enddo ! l
                  enddo ! k
            enddo ! j
      enddo ! i

      do i = 1,3
            do j = 1,3
                  Mandel(i,j) = 0.d0
                  do k = 1,3
                        do l = 1,3
      Mandel(i,j) = Mandel(i,j) + fe(k,i) * Piola(k,l) * fp(j,l)      
                        enddo ! l
                  enddo ! k
            enddo ! j
      enddo ! i

      do i = 1,3
            do j = 1,3
                  do k = 1,3
                        do l = 1,3
                              dPioladfp(i,j,k,l) = 0.d0
                              dPdf(i,j,k,l) = 0.d0
                              do m = 1,3
                                    do n = 1,3
      dPioladfp(i,j,k,l) = dPioladfp(i,j,k,l) + dpsidfe(m,n) * 
     & dfedfdfp(m,n,i,j,k,l)      
                                          do o = 1,3
                                                do p = 1,3
      dPioladfp(i,j,k,l) = dPioladfp(i,j,k,l) + dpsidfefe(m,n,o,p) * 
     & dfedf(m,n,i,j) * dfedfp(o,p,k,l)

      dPdf(i,j,k,l) = dPdf(i,j,k,l) + dpsidfefe(m,n,o,p) * 
     & dfedf(m,n,i,j) * dfedf(o,p,k,l)
                                                enddo ! p
                                          enddo ! o
                                    enddo ! n
                              enddo ! m
                        enddo ! l
                  enddo ! k
            enddo ! j
      enddo ! i

      do i = 1,3
            do j = 1,3
                  do o = 1,3
                        do p = 1,3
                              dmdfp(i,j,o,p) = 0.d0
                              dMdf(i,j,o,p) = 0.d0
                              do k = 1,3
                                    do m = 1,3
      dmdfp(i,j,o,p) = dmdfp(i,j,o,p) + dfedfp(k,i,o,p) * Piola(k,m)
     &                                      * fp(j,m)
     & + fe(k,i) * dPioladfp(k,m,o,p) * fp(j,m)
     & + fe(k,i) * Piola(k,m) * id(j,o) * id(m,p)

      dMdf(i,j,o,p) = dMdf(i,j,o,p) + dfedf(k,i,o,p) * Piola(k,m) * 
     & fp(j,m) + fe(k,i) * dPdf(k,m,o,p) * fp(j,m)
                                    enddo ! m
                              enddo ! k
                        enddo ! p
                  enddo ! o
            enddo ! j
      enddo ! i


      end subroutine
