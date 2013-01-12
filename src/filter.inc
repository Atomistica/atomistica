!
! Is atom i contained in f?
!

#define IS_EL(f, p, i)  ( iand(f, ishft(1, p%el(i))) /= 0 )
#define IS_EL2(f, e)    ( iand(f, ishft(1, e)) /= 0 )

