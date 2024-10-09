! main subroutine, which is called by python script, solves finite element calculation for one element
      subroutine getfy(d, px, py, pz, dbg, fy)

      implicit none

      real(kind=8) d(500), px, py, pz, fy

      real(kind=8) co(3,8), s(24,24), p(3,8), fext(3,8)
      real(kind=8) r(12), KK(12,12), nr, delu(12), u(3,8,2)

      integer count, dbg, i, j, info
      integer, parameter :: countmax = 50
      
Cf2py intent(in) d, px, py, pz, dbg
Cf2py intent(out) fy      

      ! p = fint
      ! s is stiffness matrix

      ! set up the problem
      co(1:2,1)   = 0.d0
      co(1,2)     = 1.d0
      co(2,2)     = 0.d0
      co(1:2,3)   = 1.d0
      co(1,4)     = 0.d0
      co(2,4)     = 1.d0
      co(3,1:4)   = 0.d0
      co(1:2,5:8) = co(1:2,1:4)
      co(3,5:8)   = 1.d0

      u = 0.d0

      ! get fext
      fext = 0.d0
      fext(3,5:8) = pz * 0.25d0
      fext(1,2:3) = px * 0.25d0
      fext(1,6:7) = px * 0.25d0
      fext(2,3:4) = py * 0.25d0
      fext(2,7:8) = py * 0.25d0

      ! do newton iteration
      nr = 1.d0 ! init just for the loop
      count = 1
      do while(nr.ge.1.d-10.and.count.le.countmax)

        ! call element routine

        call elment(d,u,co,s,p,3,3,24,.false.,fy,1.d0,8,8)

        !           dof, node
        r(1)  = fext(1,2) + p(1,2) 
        r(2)  = fext(1,3) + p(1,3) 
        r(3)  = fext(2,3) + p(2,3) 
        r(4)  = fext(2,4) + p(2,4) 
        r(5)  = fext(3,5) + p(3,5) 
        r(6)  = fext(1,6) + p(1,6) 
        r(7)  = fext(3,6) + p(3,6) 
        r(8)  = fext(1,7) + p(1,7) 
        r(9)  = fext(2,7) + p(2,7) 
        r(10) = fext(3,7) + p(3,7) 
        r(11) = fext(2,8) + p(2,8) 
        r(12) = fext(3,8) + p(3,8) 

        kk = s( (/4,7,8,11,15,16,18,19,20,21,23,24/),
     &          (/4,7,8,11,15,16,18,19,20,21,23,24/) )
            
        ! solve via Lapack
        call dposv('U', 12, 1, KK, 12, r, 12, info)

        delu = r

        ! sort back delta u for global one
        u(1,2,1) = u(1,2,1) + delu(1)
        u(1,3,1) = u(1,3,1) + delu(2)
        u(2,3,1) = u(2,3,1) + delu(3)
        u(2,4,1) = u(2,4,1) + delu(4)
        u(3,5,1) = u(3,5,1) + delu(5)
        u(1,6,1) = u(1,6,1) + delu(6)
        u(3,6,1) = u(3,6,1) + delu(7)
        u(1,7,1) = u(1,7,1) + delu(8)
        u(2,7,1) = u(2,7,1) + delu(9)
        u(3,7,1) = u(3,7,1) + delu(10)
        u(2,8,1) = u(2,8,1) + delu(11)
        u(3,8,1) = u(3,8,1) + delu(12)
        
        nr = 0.d0
        do i = 1,12
          nr = nr + r(i) * r(i)
        enddo ! i
        nr = dsqrt(nr)

        if(dbg.eq.1) then
          write(666,*) count, nr
        endif
        count = count + 1

      enddo ! while

      fy = -HUGE(fy)   ! initialise as largest possible negative number
      
      ! ask the element for fy
      call elment(d,u,co,s,p,3,3,24,.true.,fy,1.d0,8,8)

      end subroutine
