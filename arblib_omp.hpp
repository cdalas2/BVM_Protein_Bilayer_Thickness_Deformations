//
//  arblib_omp.hpp
//  
//
//  Created by CARLOS ALAS on 10/20/20.
//

#ifndef arblib_omp_hpp
#define arblib_omp_hpp

#include <arb.h>
#include <acb.h>
#include <arb_mat.h>
#include <acb_mat.h>

extern void
acb_mat_add_omp(acb_mat_t& C, const acb_mat_t& A, const acb_mat_t& B, const slong& prec, const slong& num_threads);

extern void
acb_mat_cos_entrywise_omp(acb_mat_t& B, const acb_mat_t& A, const slong& prec, const slong& num_threads);

extern void
acb_mat_div_entrywise_omp(acb_mat_t& C, const acb_mat_t& A, const acb_mat_t& B, const slong& prec, const slong& num_threads);

extern void
acb_mat_init_omp(acb_mat_t& mat, const slong& r, const slong& c, const slong& num_threads);

extern void
acb_mat_mul_entrywise_omp(acb_mat_t& C, const acb_mat_t& A, const acb_mat_t& B, const slong& prec, const slong& num_threads);

extern void
acb_mat_neg_omp(acb_mat_t& dest, const acb_mat_t& src, const slong& num_threads);

extern void
acb_mat_ones_omp(acb_mat_t& mat, const slong& num_threads);

extern void
acb_mat_scalar_div_acb_omp(acb_mat_t& B, const acb_mat_t& A, const acb_t& c, const slong& prec, const slong& num_threads);

extern void
acb_mat_scalar_div_arb_omp(acb_mat_t& B, const acb_mat_t& A, const arb_t& c, const slong& prec, const slong& num_threads);

extern void
acb_mat_scalar_div_si_omp(acb_mat_t& B, const acb_mat_t& A, const slong& c, const slong& prec, const slong& num_threads);

extern void
acb_mat_scalar_mul_acb_omp(acb_mat_t& B, const acb_mat_t& A, const acb_t& c, const slong& prec, const slong& num_threads);

extern void
acb_mat_scalar_mul_arb_omp(acb_mat_t& B, const acb_mat_t& A, const arb_t& c, const slong& prec, const slong& num_threads);

extern void
acb_mat_scalar_mul_si_omp(acb_mat_t& B, const acb_mat_t& A, const slong& c, const slong& prec, const slong& num_threads);

extern void
acb_mat_set_arb_mat_omp(acb_mat_t& dest, const arb_mat_t& src, const slong& num_threads);

extern void
acb_mat_set_omp(acb_mat_t& dest, const acb_mat_t& src, const slong& num_threads);

extern void
acb_mat_sin_cos_entrywise_omp(acb_mat_t& A, acb_mat_t& B, const acb_mat_t& C, const slong& prec, const slong& num_threads);

extern void
acb_mat_sin_entrywise_omp(acb_mat_t& B, const acb_mat_t& A, const slong& prec, const slong& num_threads);

extern void
acb_mat_sqrt_entrywise_omp(acb_mat_t& dest, const acb_mat_t& src, const slong& prec, const slong& num_threads);

extern void
acb_mat_sub_omp(acb_mat_t& C, const acb_mat_t& A, const acb_mat_t& B, const slong& prec, const slong& num_threads);

extern void
acb_mat_window_init(acb_mat_t& window, const acb_mat_t& mat, const slong& r1, const slong& c1, const slong& r2, const slong& c2);

extern void
arb_mat_add_omp(arb_mat_t& C, const arb_mat_t& A, const arb_mat_t& B, const slong& prec, const slong& num_threads);

extern void
arb_mat_cos_entrywise_omp(arb_mat_t& B, const arb_mat_t& A, const slong& prec, const slong& num_threads);

extern void arb_mat_div_entrywise_omp(arb_mat_t& C, const arb_mat_t& A, const arb_mat_t& B, const slong& prec, const slong& num_threads);

extern void
arb_mat_init_omp(arb_mat_t& mat, const slong& r, const slong& c, const slong& num_threads);

extern void
arb_mat_linspace_omp(arb_mat_t& vec, const arb_t& start_arb, const arb_t& end_arb, const slong& prec, const slong& num_threads);

extern void arb_mat_mul_entrywise_omp(arb_mat_t& C, const arb_mat_t& A, const arb_mat_t& B, const slong& prec, const slong& num_threads);

extern void
arb_mat_neg_omp(arb_mat_t& dest, const arb_mat_t& src, const slong& num_threads);

extern void
arb_mat_ones_omp(arb_mat_t& mat, const slong& num_threads);

extern void
arb_mat_scalar_addmul_si_omp(arb_mat_t& B, const arb_mat_t& A, const slong& c, const slong& prec, const slong& num_threads);

extern void
arb_mat_scalar_div_arb_omp(arb_mat_t& B, const arb_mat_t& A, const arb_t& c, const slong& prec, const slong& num_threads);

extern void
arb_mat_scalar_div_si_omp(arb_mat_t& B, const arb_mat_t& A, const slong& c, const slong& prec, const slong& num_threads);

extern void
arb_mat_scalar_mul_arb_omp(arb_mat_t& B, const arb_mat_t& A, const arb_t& c, const slong& prec, const slong& num_threads);

extern void
arb_mat_scalar_mul_si_omp(arb_mat_t& B, const arb_mat_t& A, const slong& c, const slong& prec, const slong& num_threads);

extern void
arb_mat_set_omp(arb_mat_t& dest, const arb_mat_t& src, const slong& num_threads);


extern void
arb_mat_sin_cos_entrywise_omp(arb_mat_t& A, arb_mat_t& B, const arb_mat_t& C, const slong& prec, const slong& num_threads);

extern void
arb_mat_sin_entrywise_omp(arb_mat_t& B, const arb_mat_t& A, const slong& prec, const slong& num_threads);

extern void
arb_mat_sqrt_entrywise_omp(arb_mat_t& dest, const arb_mat_t& src, const slong& prec, const slong& num_threads);

extern void
arb_mat_sub_omp(arb_mat_t& C, const arb_mat_t& A, const arb_mat_t& B, const slong& prec, const slong& num_threads);

extern void
arb_mat_window_init(arb_mat_t& window, const arb_mat_t& mat, const slong& r1, const slong& c1, const slong& r2, const slong& c2);

#endif /* arblib_omp_hpp */
