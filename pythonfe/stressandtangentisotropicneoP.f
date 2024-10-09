! calculates derivative of strain energy with respect to F^e for Neo-Hooke material
      subroutine stressandtangentisotropicneoP(f, jay, kappa, mu, 
     & Pvol, Piso, tangvol, tangiso)

      implicit none

      integer i,j,k,l

      ! in
      real*8 f(3,3), jay, kappa, mu

      ! out
      real*8 Pvol(3,3), Piso(3,3), tangvol(3,3,3,3), tangiso(3,3,3,3)

      ! calc intermediate
      real*8 jay23, finv(3,3), cbar(3,3), trcbar, dummy
      real*8 dtrcbardf(3,3), dfdfinv(3,3,3,3), lnjay, id(3,3)

      data id/ 1.d0 , 0.d0 , 0.d0 ,
     *         0.d0 , 1.d0 , 0.d0 ,
     *         0.d0 , 0.d0 , 1.d0 /

      ! preliminaries
      jay23 = jay**(-2.d0/3.d0)

      lnjay = dlog(jay)

      finv(1:3,1:3) = f(1:3,1:3)

      call invert3new(finv(1:3,1:3), dummy)

      do i = 1,3
            do j = 1,3
                  cbar(i,j) = 0.d0
                  do k = 1,3
      cbar(i,j) = cbar(i,j) + f(k,i) * f(k,j)
                  enddo ! k
            enddo ! j
      enddo ! i

      cbar = cbar * jay23

      trcbar = cbar(1,1) + cbar(2,2) + cbar(3,3)

      call dtrcbardfsub(jay23, f(1:3,1:3), finv(1:3,1:3), trcbar,
     & dtrcbardf(1:3,1:3))

      call dfinvdfsub(finv(1:3,1:3), dfdfinv(1:3,1:3,1:3,1:3))

      ! stresses and tangent
      do i = 1,3
            do j = 1,3
                  Pvol(i,j) = kappa * lnjay * finv(j,i)

                  Piso(i,j) = mu/2.d0 * dtrcbardf(i,j)
                  do k = 1,3
                        do l = 1,3
      tangvol(i,j,k,l) = kappa * (finv(l,k) * finv(j,i) + 
     & lnjay * dfdfinv(j,i,k,l))

c     is mu here correct or mu / 2?
      tangiso(i,j,k,l) = mu * (
     &      jay23 * id(i,k) * id(j,l)
     &      -f(i,j) * finv(l,k) *jay23 * 2.d0/ 3.d0
     &      - finv(j,i) * dtrcbardf(k,l) / 3.d0
     &      - trcbar/ 3.d0 * dfdfinv(j,i,k,l)
     & )
                        enddo ! l
                  enddo ! k
            enddo ! j
      enddo ! i

      end subroutine
