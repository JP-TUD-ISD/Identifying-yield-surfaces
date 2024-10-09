! calculates shapefunctions for 8 node brick element
      subroutine shpfunc(ss,det,shp,co,ndm,nel)

      implicit  none

      integer   ndm , nel , i , j , k
      real*8    rdet,det, fac1,fac2,fac3,fac4,fac5,fac6,fac7,fac8,fac9
      real*8    ss(3),shp(4,8),co(ndm,8),xs(3,3),ad(3,3)

      fac1 = 1.0d0 + ss(1)
      fac2 = 1.0d0 - ss(1)
      fac3 = 1.0d0 + ss(2)
      fac4 = 1.0d0 - ss(2)
      fac5 = 1.0d0 + ss(3)
      fac6 = 1.0d0 - ss(3)

      fac7      = 0.125d0*fac2*fac4
      fac8      = 0.125d0*fac4*fac6
      fac9      = 0.125d0*fac2*fac6
      shp(1,1)  = -fac8
      shp(1,2)  =  fac8
      shp(2,1)  = -fac9
      shp(2,4)  =  fac9
      shp(3,1)  = -fac7
      shp(3,5)  =  fac7
      shp(4,1)  =  fac7*fac6
      shp(4,5)  =  fac7*fac5

      fac7      = 0.125d0*fac1*fac3
      fac8      = 0.125d0*fac3*fac5
      fac9      = 0.125d0*fac1*fac5
      shp(1,8)  = -fac8
      shp(1,7)  =  fac8
      shp(2,6)  = -fac9
      shp(2,7)  =  fac9
      shp(3,3)  = -fac7
      shp(3,7)  =  fac7
      shp(4,3)  =  fac7*fac6
      shp(4,7)  =  fac7*fac5

      fac7      = 0.125d0*fac2*fac3
      fac8      = 0.125d0*fac4*fac5
      fac9      = 0.125d0*fac2*fac5
      shp(1,5)  = -fac8
      shp(1,6)  =  fac8
      shp(2,5)  = -fac9
      shp(2,8)  =  fac9
      shp(3,4)  = -fac7
      shp(3,8)  =  fac7
      shp(4,4)  =  fac7*fac6
      shp(4,8)  =  fac7*fac5

      fac7      = 0.125d0*fac1*fac4
      fac8      = 0.125d0*fac3*fac6
      fac9      = 0.125d0*fac1*fac6
      shp(1,4)  = -fac8
      shp(1,3)  =  fac8
      shp(2,2)  = -fac9
      shp(2,3)  =  fac9
      shp(3,2)  = -fac7
      shp(3,6)  =  fac7
      shp(4,2)  =  fac7*fac6
      shp(4,6)  =  fac7*fac5
        

      do i = 1,3
        do j = 1,3
          xs(j,i) = 0.0d0
          do k = 1,nel
            xs(j,i) = xs(j,i) + co(j,k)*shp(i,k)
          enddo ! k
        enddo ! j
      enddo ! i
      
      ad(1,1) = xs(2,2)*xs(3,3) - xs(2,3)*xs(3,2)
      ad(1,2) = xs(3,2)*xs(1,3) - xs(3,3)*xs(1,2)
      ad(1,3) = xs(1,2)*xs(2,3) - xs(1,3)*xs(2,2)
      ad(2,1) = xs(2,3)*xs(3,1) - xs(2,1)*xs(3,3)
      ad(2,2) = xs(3,3)*xs(1,1) - xs(3,1)*xs(1,3)
      ad(2,3) = xs(1,3)*xs(2,1) - xs(1,1)*xs(2,3)
      ad(3,1) = xs(2,1)*xs(3,2) - xs(2,2)*xs(3,1)
      ad(3,2) = xs(3,1)*xs(1,2) - xs(3,2)*xs(1,1)
      ad(3,3) = xs(1,1)*xs(2,2) - xs(1,2)*xs(2,1)
      det  = xs(1,1)*ad(1,1) + xs(1,2)*ad(2,1) + xs(1,3)*ad(3,1)
      rdet = 1.d0/det

      do j = 1,3
        do i = 1,3
          xs(i,j) = ad(i,j)*rdet
        enddo ! i
      enddo ! j

      do k = 1,nel
        fac7 = shp(1,k)*xs(1,1) + shp(2,k)*xs(2,1) + shp(3,k)*xs(3,1)
        fac8 = shp(1,k)*xs(1,2) + shp(2,k)*xs(2,2) + shp(3,k)*xs(3,2)
        fac9 = shp(1,k)*xs(1,3) + shp(2,k)*xs(2,3) + shp(3,k)*xs(3,3)

        
        shp(1,k) = fac7
        shp(2,k) = fac8
        shp(3,k) = fac9
      enddo ! k


      end
