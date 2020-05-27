c   Code adopted from ORAC under GPL license,
c   http://www.chim.unifi.it/orac/

      subroutine get_fftdims(nfft1,nfft2,nfft3,
     $       nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,
     $       sizfftab,sizffwrk)
      implicit none
      integer nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,
     $       nfftable,nffwork,sizfftab,sizffwrk
      integer n,nfftmax

      nfftmax = max(nfft1,nfft2,nfft3)
      nfftdim1 = nfft1
      n = nfft1/2
      if ( nfft1 .eq. 2*n )nfftdim1 = nfft1+1
      nfftdim2 = nfft2
      n = nfft2/2
      if ( nfft2 .eq. 2*n )nfftdim2 = nfft2+1
      nfftdim3 = nfft3
      n = nfft3/2
      if ( nfft3 .eq. 2*n )nfftdim3 = nfft3+1
#ifdef SGIFFT
      nfftable = 2*(nfftdim1+nfftdim2+nfftdim3+50)
      nffwork = 0
      sizfftab = nfftable
      sizffwrk  = nffwork
#endif
#ifdef CRAY
      nfftable = 2*(nfftdim1+nfftdim2+nfftdim3+50)
      nffwork = 4*nfftdim1*nfftdim2*nfftdim3
      sizfftab = nfftable
      sizffwrk  = nffwork
#else
      nfftable = 4*nfftmax + 15
      nffwork = nfftmax
      sizfftab = 3*nfftable
      sizffwrk  = 2*nfftmax
#endif
      return
      end

      subroutine fft_setup(array,fftable,ffwork,
     $      nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,
     $      nfftable,nffwork)
      implicit none

      REAL*8 array(*),fftable(*),ffwork(*)
      integer nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
      integer nfftable,nffwork,isys(4)

      integer isign,inc1,inc2,inc3
      REAL*8 scale

#ifdef SGIFFT
      call ZFFT3DI(nfft1,nfft2,nfft3,fftable)
#endif
#ifdef CRAY
      isign = 0
      scale = 1.d0
      isys(1)=3
      isys(2)=0
      isys(3)=0
      isys(4)=0
      call CCFFT3D(isign,nfft1,nfft2,nfft3,scale,array,
     $      nfftdim1,nfftdim2,array,nfftdim1,nfftdim2,fftable,
     $      ffwork,isys)
#else
      call pubz3di(nfft1,nfft2,nfft3,fftable,nfftable)
#endif
      return
      end
c-----------------------------------------------------------
      subroutine fft_forward(array,fftable,ffwork,
     $      nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,
     $      nfftable,nffwork)
      implicit none

      COMPLEX*16 array(*),ffwork(*)
      REAL*8 fftable(*)
      integer nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3

      integer isign,inc1,inc2,inc3
      REAL*8 scale
      integer nfftable,nffwork,isys(4)

      isign = 1

#ifdef SGIFFT
      call ZFFT3D(isign,nfft1,nfft2,nfft3,array,
     $   nfftdim1,nfftdim2,fftable)
#endif
#ifdef CRAY
      scale = 1.d0
      isys(1)=3
      isys(2)=0
      isys(3)=0
      isys(4)=0
      call CCFFT3D(isign,nfft1,nfft2,nfft3,scale,array,
     $      nfftdim1,nfftdim2,array,nfftdim1,nfftdim2,fftable,
     $      ffwork,isys)
#else
      call pubz3d(isign,nfft1,nfft2,nfft3,array,
     $   nfftdim1,nfftdim2,fftable,nfftable,ffwork,nffwork)
#endif
      return
      end
c-----------------------------------------------------------
      subroutine fft_back(array,fftable,ffwork,
     $      nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,
     $      nfftable,nffwork)
      implicit none

      COMPLEX*16 array(*),ffwork(*)
      REAL*8 fftable(*)
      integer nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
      integer nfftable,nffwork,isys(4)

      integer isign,inc1,inc2,inc3
      REAL*8 scale

      isign = -1

#ifdef SGIFFT
      call ZFFT3D(isign,nfft1,nfft2,nfft3,array,
     $   nfftdim1,nfftdim2,fftable)
#endif
#ifdef CRAY
      scale = 1.d0
      isys(1)=3
      isys(2)=0
      isys(3)=0
      isys(4)=0
      call CCFFT3D(isign,nfft1,nfft2,nfft3,scale,array,
     $      nfftdim1,nfftdim2,array,nfftdim1,nfftdim2,fftable,
     $      ffwork,isys)
#else
      call pubz3d(isign,nfft1,nfft2,nfft3,array,
     $   nfftdim1,nfftdim2,fftable,nfftable,ffwork,nffwork)
#endif
      return
      end
      subroutine pubz3di(n1,n2,n3,table,ntable)
      implicit none
      integer n1,n2,n3,ntable
      REAL*8 table(ntable,3)
c ntable should be 4*max(n1,n2,n3) +15


      call cffti(n1,table(1,1))
      call cffti(n2,table(1,2))
      call cffti(n3,table(1,3))

      return
      end
*****************************************************************************
      subroutine pubz3d(isign,n1,n2,n3,w,ld1,ld2,table,ntable,
     $    work,nwork)
      implicit none

      integer n1,n2,n3,ld1,ld2,isign,ntable,nwork
      COMPLEX*16 w(ld1,ld2,n3)
      COMPLEX*16 work( nwork)
      REAL*8 table(ntable,3)

      integer i,j,k
c ntable should be 4*max(n1,n2,n3) +15
c nwork should be max(n1,n2,n3)
c
c   transform along X  first ...
c
      do 100 k = 1, n3
       do 90 j = 1, n2
        do 70 i = 1,n1
          work(i) = w(i,j,k)
70      continue
        if ( isign .eq. -1) call cfftf(n1,work,table(1,1))
        if ( isign .eq. 1) call cfftb(n1,work,table(1,1))
        do 80 i = 1,n1
          w(i,j,k) = work(i)
80      continue
90     continue
100   continue
c
c   transform along Y then ...
c
      do 200 k = 1,n3
       do 190 i = 1,n1
        do 170 j = 1,n2
          work(j) = w(i,j,k)
170     continue
        if ( isign .eq. -1) call cfftf(n2,work,table(1,2))
        if ( isign .eq. 1) call cfftb(n2,work,table(1,2))
        do 180 j = 1,n2
          w(i,j,k) = work(j)
180     continue
190    continue
200   continue
c
c   transform along Z finally ...
c
      do 300 i = 1, n1
       do 290 j = 1, n2
        do 270 k = 1,n3
          work(k) = w(i,j,k)
270     continue
        if ( isign .eq. -1) call cfftf(n3,work,table(1,3))
        if ( isign .eq. 1) call cfftb(n3,work,table(1,3))
        do 280 k = 1,n3
          w(i,j,k) = work(k)
280     continue
290    continue
300   continue

      return
      end

