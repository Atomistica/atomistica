! gfortran can't handle the ## concatenation, ifort can't handle /**/ concatenation. What a pain in the ...

#ifdef __GFORTRAN__

!>
!! Inlining spline evaluation, only works for a simple_spline_t
!<

#define SPLINE_INLINE  integer :: spl_i ; real(DP) :: spl_xf, spl_B


#define SPLINE_INLINE_DEFINE(prefix, this)  real(DP) :: prefix/**/_spl_x0, prefix/**/_spl_rdx, prefix/**/_spl_y(this%n), prefix/**/_spl_coeff1(this%n-1), prefix/**/_spl_coeff2(this%n-1), prefix/**/_spl_coeff3(this%n-1), prefix/**/_spl_dcoeff1(this%n-1), prefix/**/_spl_dcoeff2(this%n-1), prefix/**/_spl_dcoeff3(this%n-1)


#define SPLINE_INLINE_OMP  private(spl_i, spl_xf, spl_B)


#define SPLINE_INLINE_OMP_DEFINE(prefix) firstprivate(prefix/**/_spl_x0, prefix/**/_spl_rdx, prefix/**/_spl_y, prefix/**/_spl_coeff1, prefix/**/_spl_coeff2, prefix/**/_spl_coeff3, prefix/**/_spl_dcoeff1, prefix/**/_spl_dcoeff2, prefix/**/_spl_dcoeff3)


#define SPLINE_INLINE_PREPARE(prefix, this)  prefix/**/_spl_x0 = this%x0 ; prefix/**/_spl_rdx = 1.0_DP/this%dx ; prefix/**/_spl_y = this%y ; prefix/**/_spl_coeff1 = this%coeff1 ; prefix/**/_spl_coeff2 = this%coeff2 ; prefix/**/_spl_coeff3 = this%coeff3 ; prefix/**/_spl_dcoeff1 = this%dcoeff1 ; prefix/**/_spl_dcoeff2 = this%dcoeff2 ; prefix/**/_spl_dcoeff3 = this%dcoeff3


#define SPLINE_FUNC(prefix, x, f)  spl_xf = (x-prefix/**/_spl_x0)*prefix/**/_spl_rdx+1.0_DP ; spl_i = spl_xf ; spl_B = spl_xf - spl_i ; f = prefix/**/_spl_y(spl_i) + spl_B*(prefix/**/_spl_coeff1(spl_i) + spl_B*(prefix/**/_spl_coeff2(spl_i) + spl_B*prefix/**/_spl_coeff3(spl_i)))


#define SPLINE_DFUNC(prefix, x, df)  spl_xf = (x-prefix/**/_spl_x0)*prefix/**/_spl_rdx+1.0_DP ; spl_i = spl_xf ; spl_B = spl_xf - spl_i ; df  = prefix/**/_spl_dcoeff1(spl_i) + spl_B*(prefix/**/_spl_dcoeff2(spl_i) + spl_B*prefix/**/_spl_dcoeff3(spl_i))


#define SPLINE_F_AND_DF(prefix, x, f, df)  spl_xf = (x-prefix/**/_spl_x0)*prefix/**/_spl_rdx+1.0_DP ; spl_i = spl_xf ; spl_B = spl_xf - spl_i ; f = prefix/**/_spl_y(spl_i) + spl_B*(prefix/**/_spl_coeff1(spl_i) + spl_B*(prefix/**/_spl_coeff2(spl_i) + spl_B*prefix/**/_spl_coeff3(spl_i))) ; df  = prefix/**/_spl_dcoeff1(spl_i) + spl_B*(prefix/**/_spl_dcoeff2(spl_i) + spl_B*prefix/**/_spl_dcoeff3(spl_i))



!>
!! Array (vectorized) spline operations
!<
#define SPLINE_INLINE_ARRAY(n)  integer :: spl_arr_i(n) ; real(DP) :: spl_arr_xf(n), spl_arr_B(n)


#define SPLINE_INLINE_ARRAY_OMP  private(spl_arr_i, spl_arr_xf, spl_arr_B)


#define SPLINE_FUNC_ARRAY(prefix, range, x, f)  spl_arr_xf(range) = (x(range)-prefix/**/_spl_x0)*prefix/**/_spl_rdx+1.0_DP ; spl_arr_i(range) = spl_arr_xf(range) ; spl_arr_B(range) = spl_arr_xf(range) - spl_arr_i(range) ; f(range) = prefix/**/_spl_y(spl_arr_i(range)) + spl_arr_B(range)*(prefix/**/_spl_coeff1(spl_arr_i(range)) + spl_arr_B(range)*(prefix/**/_spl_coeff2(spl_arr_i(range)) + spl_arr_B(range)*prefix/**/_spl_coeff3(spl_arr_i(range))))


#define SPLINE_DFUNC_ARRAY(prefix, range, x, df)  spl_arr_xf(range) = (x(range)-prefix/**/_spl_x0)*prefix/**/_spl_rdx+1.0_DP ; spl_arr_i(range) = spl_arr_xf(range) ; spl_arr_B(range) = spl_arr_xf(range) - spl_arr_i(range) ; df(range)  = prefix/**/_spl_dcoeff1(spl_arr_i(range)) + spl_arr_B(range)*(prefix/**/_spl_dcoeff2(spl_arr_i(range)) + spl_arr_B(range)*prefix/**/_spl_dcoeff3(spl_arr_i(range)))


#define SPLINE_F_AND_DF_ARRAY(prefix, range, x, f, df)  spl_arr_xf(range) = (x(range)-prefix/**/_spl_x0)*prefix/**/_spl_rdx+1.0_DP ; spl_arr_i(range) = spl_arr_xf(range) ; spl_arr_B(range) = spl_arr_xf(range) - spl_arr_i(range) ; f(range) = prefix/**/_spl_y(spl_arr_i(range)) + spl_arr_B(range)*(prefix/**/_spl_coeff1(spl_arr_i(range)) + spl_arr_B(range)*(prefix/**/_spl_coeff2(spl_arr_i(range)) + spl_arr_B(range)*prefix/**/_spl_coeff3(spl_arr_i(range)))) ; df(range)  = prefix/**/_spl_dcoeff1(spl_arr_i(range)) + spl_arr_B(range)*(prefix/**/_spl_dcoeff2(spl_arr_i(range)) + spl_arr_B(range)*prefix/**/_spl_dcoeff3(spl_arr_i(range)))

#else

!>
!! Inlining spline evaluation, only works for a simple_spline_t
!<

#define SPLINE_INLINE  integer :: spl_i ; real(DP) :: spl_xf, spl_B


#define SPLINE_INLINE_DEFINE(prefix, this)  real(DP) :: prefix ## _spl_x0, prefix ## _spl_rdx, prefix ## _spl_y(this%n), prefix ## _spl_coeff1(this%n-1), prefix ## _spl_coeff2(this%n-1), prefix ## _spl_coeff3(this%n-1), prefix ## _spl_dcoeff1(this%n-1), prefix ## _spl_dcoeff2(this%n-1), prefix ## _spl_dcoeff3(this%n-1)


#define SPLINE_INLINE_OMP  private(spl_i, spl_xf, spl_B)


#define SPLINE_INLINE_OMP_DEFINE(prefix) firstprivate(prefix ## _spl_x0, prefix ## _spl_rdx, prefix ## _spl_y, prefix ## _spl_coeff1, prefix ## _spl_coeff2, prefix ## _spl_coeff3, prefix ## _spl_dcoeff1, prefix ## _spl_dcoeff2, prefix ## _spl_dcoeff3)


#define SPLINE_INLINE_PREPARE(prefix, this)  prefix ## _spl_x0 = this%x0 ; prefix ## _spl_rdx = 1.0_DP/this%dx ; prefix ## _spl_y = this%y ; prefix ## _spl_coeff1 = this%coeff1 ; prefix ## _spl_coeff2 = this%coeff2 ; prefix ## _spl_coeff3 = this%coeff3 ; prefix ## _spl_dcoeff1 = this%dcoeff1 ; prefix ## _spl_dcoeff2 = this%dcoeff2 ; prefix ## _spl_dcoeff3 = this%dcoeff3


#define SPLINE_FUNC(prefix, x, f)  spl_xf = (x-prefix ## _spl_x0)*prefix ## _spl_rdx+1.0_DP ; spl_i = spl_xf ; spl_B = spl_xf - spl_i ; f = prefix ## _spl_y(spl_i) + spl_B*(prefix ## _spl_coeff1(spl_i) + spl_B*(prefix ## _spl_coeff2(spl_i) + spl_B*prefix ## _spl_coeff3(spl_i)))


#define SPLINE_DFUNC(prefix, x, df)  spl_xf = (x-prefix ## _spl_x0)*prefix ## _spl_rdx+1.0_DP ; spl_i = spl_xf ; spl_B = spl_xf - spl_i ; df  = prefix ## _spl_dcoeff1(spl_i) + spl_B*(prefix ## _spl_dcoeff2(spl_i) + spl_B*prefix ## _spl_dcoeff3(spl_i))


#define SPLINE_F_AND_DF(prefix, x, f, df)  spl_xf = (x-prefix ## _spl_x0)*prefix ## _spl_rdx+1.0_DP ; spl_i = spl_xf ; spl_B = spl_xf - spl_i ; f = prefix ## _spl_y(spl_i) + spl_B*(prefix ## _spl_coeff1(spl_i) + spl_B*(prefix ## _spl_coeff2(spl_i) + spl_B*prefix ## _spl_coeff3(spl_i))) ; df  = prefix ## _spl_dcoeff1(spl_i) + spl_B*(prefix ## _spl_dcoeff2(spl_i) + spl_B*prefix ## _spl_dcoeff3(spl_i))



!>
!! Array (vectorized) spline operations
!<
#define SPLINE_INLINE_ARRAY(n)  integer :: spl_arr_i(n) ; real(DP) :: spl_arr_xf(n), spl_arr_B(n)


#define SPLINE_INLINE_ARRAY_OMP  private(spl_arr_i, spl_arr_xf, spl_arr_B)


#define SPLINE_FUNC_ARRAY(prefix, range, x, f)  spl_arr_xf(range) = (x(range)-prefix ## _spl_x0)*prefix ## _spl_rdx+1.0_DP ; spl_arr_i(range) = spl_arr_xf(range) ; spl_arr_B(range) = spl_arr_xf(range) - spl_arr_i(range) ; f(range) = prefix ## _spl_y(spl_arr_i(range)) + spl_arr_B(range)*(prefix ## _spl_coeff1(spl_arr_i(range)) + spl_arr_B(range)*(prefix ## _spl_coeff2(spl_arr_i(range)) + spl_arr_B(range)*prefix ## _spl_coeff3(spl_arr_i(range))))


#define SPLINE_DFUNC_ARRAY(prefix, range, x, df)  spl_arr_xf(range) = (x(range)-prefix ## _spl_x0)*prefix ## _spl_rdx+1.0_DP ; spl_arr_i(range) = spl_arr_xf(range) ; spl_arr_B(range) = spl_arr_xf(range) - spl_arr_i(range) ; df(range)  = prefix ## _spl_dcoeff1(spl_arr_i(range)) + spl_arr_B(range)*(prefix ## _spl_dcoeff2(spl_arr_i(range)) + spl_arr_B(range)*prefix ## _spl_dcoeff3(spl_arr_i(range)))


#define SPLINE_F_AND_DF_ARRAY(prefix, range, x, f, df)  spl_arr_xf(range) = (x(range)-prefix ## _spl_x0)*prefix ## _spl_rdx+1.0_DP ; spl_arr_i(range) = spl_arr_xf(range) ; spl_arr_B(range) = spl_arr_xf(range) - spl_arr_i(range) ; f(range) = prefix ## _spl_y(spl_arr_i(range)) + spl_arr_B(range)*(prefix ## _spl_coeff1(spl_arr_i(range)) + spl_arr_B(range)*(prefix ## _spl_coeff2(spl_arr_i(range)) + spl_arr_B(range)*prefix ## _spl_coeff3(spl_arr_i(range)))) ; df(range)  = prefix ## _spl_dcoeff1(spl_arr_i(range)) + spl_arr_B(range)*(prefix ## _spl_dcoeff2(spl_arr_i(range)) + spl_arr_B(range)*prefix ## _spl_dcoeff3(spl_arr_i(range)))

#endif