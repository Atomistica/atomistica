C R250.F77     The R250 Pseudo-random number generator
C
C algorithm from:
C Kirkpatrick, S., and E. Stoll, 1981; A Very Fast Shift-Register
C Sequence Random Number Generator, Journal of Computational Physics,
C V. 40. p. 517
C 
C see also:
C Maier, W.L., 1991; A Fast Pseudo Random Number Generator,
C                    Dr. Dobb's Journal, May, pp. 152 - 157
C
C 
C Uses the Linear Congruential Method,
C the "minimal standard generator"
C Park & Miller, 1988, Comm of the ACM, 31(10), pp. 1192-1201
C for initialization
C
C
C For a review of BOTH of these generators, see:
C Carter, E.F, 1994; Generation and Application of Random Numbers,
C Forth Dimensions, Vol. XVI, Numbers 1,2 May/June, July/August
C
C
C $Author:   skip  $
C $Workfile:   r250.f  $
C $Revision:   1.1  $
C $Date:   07 Nov 1996 01:23:06  $
C
C Modified: 12 Jan 2009 for OpenMP, i.e. RNG for each thread.
C
C ===================================================================
C
      Function r250_lcmrand(ix, x)
C     The minimal standard PRNG for 31 bit unsigned integers
C     designed with automatic overflow protection  
C     uses ix as the seed value if it is greater than zero
C     otherwise it is ignored
      Integer*4 ix
      Integer*4 a, b, m, q, r
      Integer*4 hi, lo, test
      Integer*4 x
      Parameter (a = 16807, b = 0, m = 2147483647)
      Parameter (q = 127773, r = 2836)
C
      If ( ix .gt. 0 ) x = ix
      
      hi = x / q
      lo = mod( x, q )
      test = a * lo - r * hi
      if ( test .gt. 0 ) then
          x = test
      else
          x = test + m
      endif
      
      r250_lcmrand = x
      return
      End
  

C ===================================================================
C
C  R250, call R250Init with the desired initial seed BEFORE
C  the first invocation of IRAND()
C
C ===================================================================

      Subroutine r250_init(iseed,indexf,indexb,buffer)
      Integer*4 k, mask, msb
      Integer*4 indexf, indexb, buffer(250)
      Integer*4 ms_bit, all_bits, half_range, step
      Integer*4 x
      Parameter ( ms_bit = Z'40000000')
      Parameter ( half_range = Z'20000000' )
      Parameter ( all_bits = Z'7FFFFFFF' )
      Parameter ( step = 7 )
C
      indexf = 1
      indexb = 104
      k = iseed
      Do 10 i = 1, 250
	  buffer(i) = r250_lcmrand( k, x )
	  k = -1
  10  EndDo
      Do 20 i = 1, 250
         if ( r250_lcmrand( -1, x ) .gt. half_range ) then
	     buffer(i) = ior( buffer(i), ms_bit )
	 endif
 20   EndDo

      msb = ms_bit
      mask = all_bits

      Do 30 i = 0,30
        k = step * i + 4
        buffer(k) = iand( buffer(k), mask )
        buffer(k) = ior( buffer(k), msb )
        msb = msb / 2
        mask = mask / 2
  30  EndDo
      
      Return
      END



      integer*4 Function r250_irand(indexf,indexb,buffer)
C     R250 PRNG, run after R250_Init
      Integer*4 newrand
      Integer*4 indexf, indexb, buffer(250)

      newrand = ieor( buffer(indexf), buffer(indexb) )
      buffer(indexf) = newrand

      indexf = indexf + 1
      if ( indexf .gt. 250 ) indexf = 1

      indexb = indexb + 1
      if ( indexb .gt. 250 ) indexb = 1
	  

      r250_irand = newrand
      return
      End
  
