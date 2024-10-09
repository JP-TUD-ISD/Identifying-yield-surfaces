! finite element formulation, balance of linear momentum expressed in first Piola-Kirchhoff stresses
      subroutine elment(d,ul,xl,s,p,ndf,ndm,nst,fyflag,fy,dt,nen,nel)

      implicit   none

      integer, parameter :: maxguass = 8
      integer, parameter :: maxnel = 10
      integer ndf , ndm , nst 
      logical fyflag
      real*8  d(500)  , ul(ndf,nen,2)  , xl(ndm,nen)
      real*8  p(ndf,nen) , s(nst,nst)  
      
      ! integer 
      integer i,j,k,l,ii,jj,kk,ll,mm,nn,ngauss,ijump,jjump,point
      
      ! switches
      integer msw,fibsw
      
      ! output
      real*8 quant(50,maxguass)
      
      ! history
      integer hist
      real*8, allocatable :: hn(:),h1(:)
      
      ! variables from finite Element Method
      real*8 shp(4,maxnel,maxguass),f(3,3,2,maxguass),finv(3,3,maxguass)
      real*8 detf(2,maxguass),dvol0(maxguass)
      real*8 dvol(maxguass),remshp(maxguass),df(3,3,maxguass)
      real*8 resi(ndf,maxnel),Jay(maxguass)
      ! remark Jay means determinant of Jacobi matrix of the ansatz functions, wrt to local coordinates
      real*8 sg(4,maxguass),refshp(4,maxnel,maxguass),dt2,sv(5,maxguass)

      ! quantities at local gaussian point, coming from material model and material library
      real*8 Piolal(3,3),tangPiolal(3,3,3,3),rho0l,rhol, dplot

      ! quantities at element level from gaussian points
      real*8 Piola(3,3,maxguass),tangPiola(3,3,3,3,maxguass)
      real*8 rho0(maxguass),rho(maxguass)
      real*8 sigplot(6),sigsig(6,maxguass)

      real*8 udotdot(3,maxguass),udot(3,maxguass)

      ! constants
      real*8 bodyforce(3),id(3,3),pi,hnd(1),h1d(1)
      integer ij(3,3),kl(3,3)

      ! for testing only
      real*8 sigtest(6),ctest(6,6)

      ! for viscosity
      integer nbranches, hbr, br

      ! for transitioning
      real(kind=8) dt, fy
      integer nel,nen

      ! internal 
      integer*8 xpoint, xp2
      integer xlen,xprec, xl2, xprec2
      logical xflag,ualloc,setvar, xfl2
      
      ! area for data statements
      data ngauss / 8 /
      !data bodyforce / 0.d0 , 0.d0 , -9.81d0 / ! -9.81d0
      data ij/ 1 , 4 , 6 ,
     *         4 , 2 , 5 ,
     *         6 , 5 , 3 /  
      data kl/ 1 , 7 , 9 ,
     *         4 , 2 , 8 ,
     *         6 , 5 , 3 / 
      data id/ 1.d0 , 0.d0 , 0.d0 ,
     *         0.d0 , 1.d0 , 0.d0 ,
     *         0.d0 , 0.d0 , 1.d0 /
      data pi/3.14159265359d0/

      dt2 = dt
      if(dt2.le.1.d-10) dt2 = 1.d-10
      
      msw             = int(d(1))
      hist            = int(d(3)) 
      bodyforce(1:3)  = d(4:6)
      rho0l           = d(2)

      ! start with zero the FEM arrays for global solver 
      p(1:ndf,1:nen) = 0.d0
      s(1:nst,1:nst) = 0.d0
      
      call sweight(ngauss,sg) 

      ! calculate the shape functions
      do point = 1,ngauss
              
      call shpfunc(sg(1:3,point),Jay(point),shp(1:4,1:8,point),
     2 xl(1:ndm,1:8),ndm,nel)
       dvol0(point) = sg(4,point) * Jay(point) ! this is gaussian weight * det(J^e)      
      enddo ! point

              do point = 1,ngauss
                  
                  call calcfs(shp(1:4,1:nel,point), 
     2            ul(1:ndf,1:nen,1:2), 
     2 f(1:3,1:3,1:2,point), finv(1:3,1:3,point),df(1:3,1:3,point)
     2 ,detf(1:2,point),
     3 ndf,nel,nen)
                  
               ! compute the current volume
                  
               dvol(point) = dvol0(point) * detf(1,point)
               
               ! compute shape functions at the current configuration
               do k = 1, nel
                   do j = 1,3
                       remshp(j) = 0.d0
                       do i = 1,3 ! multiply the shape functions in ref config with F^inv to transform the gradient into current config and store in remshp
              remshp(j) = remshp(j) + finv(i,j,point) * shp(i,k,point)
                       enddo ! i
                   enddo ! j
                   do j = 1,3 ! replace the gradient with basis in ref config with gradient in current config
                       refshp(j,k,point) = shp(j,k,point) ! save the reference configuration shape functions
                       shp(j,k,point) = remshp(j)
                   enddo ! j
                   refshp(4,k,point) = shp(4,k,point) ! 4th entry should not be changed ... readability improvement
               enddo ! k 
                        
      udot(1:3,point) = 0.d0
              
      enddo ! point
              
          
      ! now all quantities from the global solution are prepared to be plugged into the next stiffness and residual
      do point = 1,ngauss
          
            
      quant(:,point) = 0.d0
      if(hist.gt.0)then
      
      ! here call the new matlib
      call Piolamatlib(d,f(1:3,1:3,1:2,point), detf(1:2,point), 
     2 finv(1:3,1:3,point), Piolal(1:3,1:3), tangPiolal(1:3,1:3,1:3,1:3)
     2 , hist, hn(1:hist), h1(1:hist), fyflag, msw, point, dplot, fy)


      else
      ! here also call the new matlib but with different h1 and hn
      call Piolamatlib(d,f(1:3,1:3,1:2,point), detf(1:2,point), 
     2 finv(1:3,1:3,point), Piolal(1:3,1:3), tangPiolal(1:3,1:3,1:3,1:3)
     2 , 1, hnd(1:1), h1d(1:1), fyflag, msw, point, dplot, fy)

      endif
      
      
      
! here the integration part can already be done
      Piola(1:3,1:3,point) = Piolal(1:3,1:3) * dvol0(point)
      tangPiola(1:3,1:3,1:3,1:3,point) = 
     3 tangPiolal(1:3,1:3,1:3,1:3) * dvol0(point)
      rho0(point) = rho0l 
      rhol = rho0l / detf(1,point)
      rho(point) = rhol 

      
      enddo ! point


          do point = 1,ngauss
          do i = 1,nel
              resi(1:ndf,i) = 0.d0 ! zero for sums
          ! here residual 1 is calculated, which is used for node displacements
                  do ii = 1,3
                      resi(ii,i) = resi(ii,i) 
     2 - refshp(4,i,point) * rho0(point) * bodyforce(ii) * dvol0(point) ! ref volume because different form of balance of linear momentum
                      
                      do kk = 1,3
                          
                          resi(ii,i) = resi(ii,i) + 
     2 refshp(kk,i,point) * Piola(ii,kk,point) ! Piola already integrated over reference volume
                          
                      enddo ! kk
                  enddo ! ii
                  
          p(1:ndf,i) = p(1:ndf,i) - resi(1:ndf,i)
          
          enddo ! i
          enddo ! point
          
          
              do point = 1,ngauss
                  
                  ijump = 0
                  do i = 1,nel
                      jjump = 0
                      do j = 1,nel
                          do ii = 1,3
                              do jj = 1,3                                  
                                  do kk = 1,3
                                      do ll = 1,3
                                          
      ! this part below is the part from the linearization of the kirchoff stresses
      s(ijump+ii,jjump+jj) = s(ijump+ii,jjump+jj) 
     2 + refshp(kk,i,point) * tangPiola(ii,kk,jj,ll,point) 
     2 * refshp(ll,j,point)
                                      enddo ! ll
                                  enddo ! kk
                              enddo ! jj
                          enddo ! ii
                          jjump = jjump + ndf
                      enddo ! j
                      ijump = ijump + ndf
                  enddo ! i
              enddo ! point
              

         
      if(allocated(hn)) deallocate(hn)
      if(allocated(h1)) deallocate(h1)
          


      end subroutine elment
