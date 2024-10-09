! calculates rotated vectors A and B, alpha contains angles for rotation around x_1, x_2 and x_3 axis
      subroutine rotateAB(alpha,A,B)
      implicit none
      real*8 alpha(3),A(3),B(3),alpx,alpy,alpz
      real*8 rotx(3,3),roty(3,3),rotz(3,3),rot(3,3)
      real*8 an(3),bn(3)
      alpx=alpha(1)
      alpy=alpha(2)
      alpz=alpha(3)
      
      rotx(1:3,1:3) = 0.d0
      rotx(1,1)     = 1.d0
      rotx(2,2)     = dcos(alpx)
      rotx(3,3)     = dcos(alpx)
      rotx(2,3)     = -dsin(alpx)
      rotx(3,2)     = dsin(alpx)
      
      roty(1:3,1:3) = 0.d0
      roty(2,2)     = 1.d0
      roty(1,1)     = dcos(alpy)
      roty(3,3)     = dcos(alpy)
      roty(3,1)     = -dsin(alpy)
      roty(1,3)     = dsin(alpy)
      
      rotz(1:3,1:3) = 0.d0
      rotz(3,3)     = 1.d0
      rotz(2,2)     = dcos(alpz)
      rotz(1,1)     = dcos(alpz)
      rotz(1,2)     = -dsin(alpz)
      rotz(2,1)     = dsin(alpz)
      
      rot = matmul(rotz,matmul(roty,rotx))
      an  = a
      bn  = b
      a   = matmul(rot,an)
      b   = matmul(rot,bn)
      end subroutine    
