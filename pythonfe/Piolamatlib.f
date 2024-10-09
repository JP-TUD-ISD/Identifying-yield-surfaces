! material library called from element, can be extended to include previous history
      subroutine Piolamatlib(d,f,detf,finv,Piola,tangPiola,nh,hn,h1,
     2 fyflag, msw,gp, dplot, fyo)

      implicit none

      integer nh,msw,mls,gp,i,j,k,l,kl(3,3)
      logical fyflag
      real*8 d(*),f(3,3,2),detf(2),finv(3,3),Piola(3,3)
      real*8 tangPiola(3,3,3,3),hn(nh),h1(nh)
      real*8 sig(9),tang(9,9),zeros(6),zeros_m(3,3), dplot, fyo

      data kl/ 1 , 7 , 9 ,
     *         4 , 2 , 8 ,
     *         6 , 5 , 3 / 

      if(msw.eq.6969) then

      call plastwo(d, f, detf, fyflag, Piola, tangPiola, fyo)

      endif

      end
