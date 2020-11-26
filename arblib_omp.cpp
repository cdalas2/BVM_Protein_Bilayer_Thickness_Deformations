//
//  arblib_omp.cpp
//  
//
//  Created by CARLOS ALAS on 10/20/20.
//

#include "arblib_omp.hpp"
#include <arb_mat.h>
#include <arb.h>
#include <acb_mat.h>
#include <acb.h>
#include <omp.h>

void
acb_mat_add_omp(acb_mat_t& C, const acb_mat_t& A, const acb_mat_t& B, const slong& prec, const slong& num_threads)
{
    if ( (acb_mat_nrows(B) != acb_mat_nrows(A)) ||
         (acb_mat_ncols(B) != acb_mat_ncols(A)) )
    {
        flint_printf("acb_mat_add_omp: incompatible dimensions\n");
        flint_abort();
    }
    
        
    slong r, c;

    r = acb_mat_nrows(C);
    c = acb_mat_ncols(C);

    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < r; i++)
        for (slong j = 0; j < c; j++)
        acb_add(acb_mat_entry(C,i,j),acb_mat_entry(A,i,j),acb_mat_entry(B,i,j),prec);
        
        flint_cleanup();
    }
};


void
acb_mat_cos_entrywise_omp(acb_mat_t& B, const acb_mat_t& A, const slong& prec, const slong& num_threads)
{
    if ( (acb_mat_nrows(A) != acb_mat_nrows(B)) ||
         (acb_mat_ncols(A) != acb_mat_ncols(B)) )
    {
        flint_printf("acb_mat_cos_entrywise_omp: incompatible dimensions\n");
        flint_abort();
    }
    
        
    slong r, c;

    r = acb_mat_nrows(B);
    c = acb_mat_ncols(B);

    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < r; i++)
        for (slong j = 0; j < c; j++)
        acb_cos(acb_mat_entry(B,i,j),acb_mat_entry(A,i,j),prec);
        
        flint_cleanup();
    }
};


void
acb_mat_div_entrywise_omp(acb_mat_t& C, const acb_mat_t& A, const acb_mat_t& B, const slong& prec, const slong& num_threads)
{

    if (acb_mat_nrows(A) != acb_mat_nrows(B) ||
        acb_mat_ncols(A) != acb_mat_ncols(B))
    {
        flint_printf("acb_mat_div_entrywise_omp: incompatible dimensions\n");
        flint_abort();
    }
    
    slong r = acb_mat_nrows(A);
    slong c = acb_mat_ncols(A);
    
    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < r; i++)
        {
            for (slong j = 0; j < c; j++)
            {
                acb_div(arb_mat_entry(C, i, j),
                        acb_mat_entry(A, i, j),
                        acb_mat_entry(B, i, j), prec);
            }
        }
        flint_cleanup();
    }
};


void
acb_mat_init_omp(acb_mat_t& mat, const slong& r, const slong& c, const slong& num_threads)
{
    if (r != 0 && c != 0)
    {
        mat->entries = _acb_vec_init(r * c);
        mat->rows = (acb_ptr *) flint_malloc(r * sizeof(acb_ptr));
        omp_set_num_threads(num_threads);
        #pragma omp parallel
        {
            #pragma omp for
            for (slong i = 0; i < r; i++)
            mat->rows[i] = mat->entries + i * c;
            
            flint_cleanup();
        }
    }
    else
        mat->entries = NULL;

    mat->r = r;
    mat->c = c;
};

void acb_mat_mul_entrywise_omp(acb_mat_t& C, const acb_mat_t& A, const acb_mat_t& B, const slong& prec, const slong& num_threads)
{
    slong i, j;

    if (acb_mat_nrows(A) != acb_mat_nrows(B) ||
        acb_mat_ncols(A) != acb_mat_ncols(B))
    {
        flint_printf("acb_mat_mul_entrywise_omp: incompatible dimensions\n");
        flint_abort();
    }
    
    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (i = 0; i < acb_mat_nrows(A); i++)
        {
            for (j = 0; j < acb_mat_ncols(A); j++)
            {
                acb_mul(acb_mat_entry(C, i, j),
                        acb_mat_entry(A, i, j),
                        acb_mat_entry(B, i, j), prec);
            }
        }
        flint_cleanup();
    }
};


void
acb_mat_neg_omp(acb_mat_t& dest, const acb_mat_t& src, const slong& num_threads)
{
    if ( (acb_mat_nrows(dest) != acb_mat_nrows(src)) ||
         (acb_mat_ncols(dest) != acb_mat_ncols(src)) )
    {
        flint_printf("acb_mat_neg_omp: incompatible dimensions\n");
        flint_abort();
    }
    
    slong R, C;

    R = acb_mat_nrows(src);
    C = acb_mat_ncols(src);

    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < R; i++)
        for (slong j = 0; j < C; j++)
        acb_neg(acb_mat_entry(dest, i, j), acb_mat_entry(src, i, j));
        
        flint_cleanup();
    }
};


void
acb_mat_ones_omp(acb_mat_t& mat, const slong& num_threads)
{
    slong R, C;

    R = acb_mat_nrows(mat);
    C = acb_mat_ncols(mat);

    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < R; i++)
        for (slong j = 0; j < C; j++)
        acb_one(acb_mat_entry(mat, i, j));
        
        flint_cleanup();
    }
};


void
acb_mat_scalar_div_acb_omp(acb_mat_t& B, const acb_mat_t& A, const acb_t& c, const slong& prec, const slong& num_threads)
{
    if ( (acb_mat_nrows(B) != acb_mat_nrows(A)) ||
         (acb_mat_ncols(B) != acb_mat_ncols(A)) )
    {
        flint_printf("acb_mat_scalar_div_acb_omp: incompatible dimensions\n");
        flint_abort();
    }
    
        
    slong R, C;

    R = acb_mat_nrows(B);
    C = acb_mat_ncols(B);

    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < R; i++)
        for (slong j = 0; j < C; j++)
        acb_div(acb_mat_entry(B,i,j),acb_mat_entry(A,i,j),c,prec);
        
        flint_cleanup();
    }
};


void
acb_mat_scalar_div_arb_omp(acb_mat_t& B, const acb_mat_t& A, const arb_t& c, const slong& prec, const slong& num_threads)
{
    if ( (acb_mat_nrows(B) != acb_mat_nrows(A)) ||
         (acb_mat_ncols(B) != acb_mat_ncols(A)) )
    {
        flint_printf("acb_mat_scalar_div_arb_omp: incompatible dimensions\n");
        flint_abort();
    }
    
        
    slong R, C;

    R = acb_mat_nrows(B);
    C = acb_mat_ncols(B);

    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < R; i++)
        for (slong j = 0; j < C; j++)
        acb_div_arb(acb_mat_entry(B,i,j),acb_mat_entry(A,i,j),c,prec);
        
        flint_cleanup();
    }
};


void
acb_mat_scalar_div_si_omp(acb_mat_t& B, const acb_mat_t& A, const slong& c, const slong& prec, const slong& num_threads)
{
    if ( (acb_mat_nrows(B) != acb_mat_nrows(A)) ||
         (acb_mat_ncols(B) != acb_mat_ncols(A)) )
    {
        flint_printf("acb_mat_scalar_div_si_omp: incompatible dimensions\n");
        flint_abort();
    }
    
        
    slong R, C;

    R = acb_mat_nrows(A);
    C = acb_mat_ncols(A);

    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < R; i++)
        for (slong j = 0; j < C; j++)
        acb_div_si(acb_mat_entry(B,i,j),acb_mat_entry(A,i,j),c,prec);
        
        flint_cleanup();
    }
};


void
acb_mat_scalar_mul_acb_omp(acb_mat_t& B, const acb_mat_t& A, const acb_t& c, const slong& prec, const slong& num_threads)
{
    if ( (acb_mat_nrows(B) != acb_mat_nrows(A)) ||
         (acb_mat_ncols(B) != acb_mat_ncols(A)) )
    {
        flint_printf("acb_mat_scalar_mul_acb_omp: incompatible dimensions\n");
        flint_abort();
    }
    
        
    slong R, C;

    R = acb_mat_nrows(B);
    C = acb_mat_ncols(B);

    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < R; i++)
        for (slong j = 0; j < C; j++)
        acb_mul(acb_mat_entry(B,i,j),acb_mat_entry(A,i,j),c,prec);
        
        flint_cleanup();
    }
};


void
acb_mat_scalar_mul_arb_omp(acb_mat_t& B, const acb_mat_t& A, const arb_t& c, const slong& prec, const slong& num_threads)
{
    if ( (acb_mat_nrows(B) != acb_mat_nrows(A)) ||
         (acb_mat_ncols(B) != acb_mat_ncols(A)) )
    {
        flint_printf("acb_mat_scalar_mul_arb_omp: incompatible dimensions\n");
        flint_abort();
    }
    
        
    slong R, C;

    R = acb_mat_nrows(B);
    C = acb_mat_ncols(B);

    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < R; i++)
        for (slong j = 0; j < C; j++)
        acb_mul_arb(acb_mat_entry(B,i,j),acb_mat_entry(A,i,j),c,prec);
        
        flint_cleanup();
    }
};


void
acb_mat_scalar_mul_si_omp(acb_mat_t& B, const acb_mat_t& A, const slong& c, const slong& prec, const slong& num_threads)
{
    if ( (acb_mat_nrows(B) != acb_mat_nrows(A)) ||
         (acb_mat_ncols(B) != acb_mat_ncols(A)) )
    {
        flint_printf("acb_mat_scalar_mul_si_omp: incompatible dimensions\n");
        flint_abort();
    }
    
        
    slong R, C;

    R = acb_mat_nrows(A);
    C = acb_mat_ncols(A);

    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < R; i++)
        for (slong j = 0; j < C; j++)
        acb_mul_si(acb_mat_entry(B,i,j),acb_mat_entry(A,i,j),c,prec);
        
        flint_cleanup();
    }
};


void
acb_mat_set_arb_mat_omp(acb_mat_t& dest, const arb_mat_t& src, const slong& num_threads)
{
    slong R, C;

    R = acb_mat_nrows(dest);
    C = acb_mat_ncols(dest);

    omp_set_num_threads(num_threads);
    
    if (arb_mat_ncols(dest) != 0)
    {
        #pragma omp parallel
        {
            #pragma omp for
            for (slong i = 0; i < R; i++)
            for (slong j = 0; j < C; j++)
            acb_set_arb(acb_mat_entry(dest, i, j), arb_mat_entry(src, i, j));
            
            flint_cleanup();
        }
    }
};


void
acb_mat_set_omp(acb_mat_t& dest, const acb_mat_t& src, const slong& num_threads)
{
    slong R, C;

    R = acb_mat_nrows(dest);
    C = acb_mat_ncols(dest);

    omp_set_num_threads(num_threads);
    
    if (dest != src && acb_mat_ncols(src) != 0)
    {
        #pragma omp parallel
        {
            #pragma omp for
            for (slong i = 0; i < R; i++)
            for (slong j = 0; j < C; j++)
            acb_set(acb_mat_entry(dest, i, j), acb_mat_entry(src, i, j));
            
            flint_cleanup();
        }
    }
};


void
acb_mat_sin_cos_entrywise_omp(acb_mat_t& A, acb_mat_t& B, const acb_mat_t& C, const slong& prec, const slong& num_threads)
{
    if ( (acb_mat_nrows(C) != acb_mat_nrows(A)) ||
         (acb_mat_ncols(C) != acb_mat_ncols(A)) ||
         (acb_mat_nrows(C) != acb_mat_nrows(B)) ||
         (acb_mat_ncols(C) != acb_mat_ncols(B)) )
    {
        flint_printf("acb_mat_sin_cos_entrywise_omp: incompatible dimensions\n");
        flint_abort();
    }
    
        
    slong r, c;

    r = acb_mat_nrows(C);
    c = acb_mat_ncols(C);

    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < r; i++)
        for (slong j = 0; j < c; j++)
        acb_sin_cos(acb_mat_entry(A,i,j),acb_mat_entry(B,i,j),acb_mat_entry(C,i,j),prec);
        
        flint_cleanup();
    }
};


void
acb_mat_sin_entrywise_omp(acb_mat_t& B, const acb_mat_t& A, const slong& prec, const slong& num_threads)
{
    if ( (acb_mat_nrows(A) != acb_mat_nrows(B)) ||
         (acb_mat_ncols(A) != acb_mat_ncols(B)) )
    {
        flint_printf("acb_mat_sin_entrywise_omp: incompatible dimensions\n");
        flint_abort();
    }
    
        
    slong r, c;

    r = acb_mat_nrows(B);
    c = acb_mat_ncols(B);

    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < r; i++)
        for (slong j = 0; j < c; j++)
        acb_sin(acb_mat_entry(B,i,j),acb_mat_entry(A,i,j),prec);
        
        flint_cleanup();
    }
};


void
acb_mat_sqrt_entrywise_omp(acb_mat_t& dest, const acb_mat_t& src, const slong& prec, const slong& num_threads)
{
    if ( (acb_mat_nrows(dest) != acb_mat_nrows(src)) ||
         (acb_mat_ncols(dest) != acb_mat_ncols(src)) )
    {
        flint_printf("acb_mat_sqrt_entrywise_omp: incompatible dimensions\n");
        flint_abort();
    }
    
    slong R, C;

    R = acb_mat_nrows(src);
    C = acb_mat_ncols(src);

    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < R; i++)
        for (slong j = 0; j < C; j++)
        acb_sqrt(acb_mat_entry(dest, i, j), acb_mat_entry(src, i, j),prec);
        
        flint_cleanup();
    }
};


void
acb_mat_sub_omp(acb_mat_t& C, const acb_mat_t& A, const acb_mat_t& B, const slong& prec, const slong& num_threads)
{
    if ( (acb_mat_nrows(B) != acb_mat_nrows(A)) ||
         (acb_mat_ncols(B) != acb_mat_ncols(A)) )
    {
        flint_printf("acb_mat_sub_omp: incompatible dimensions\n");
        flint_abort();
    }
    
        
    slong r, c;

    r = acb_mat_nrows(C);
    c = acb_mat_ncols(C);

    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < r; i++)
        for (slong j = 0; j < c; j++)
        acb_sub(acb_mat_entry(C,i,j),acb_mat_entry(A,i,j),acb_mat_entry(B,i,j),prec);
        
        flint_cleanup();
    }
};


void
acb_mat_window_init(acb_mat_t& window, const acb_mat_t& mat, const slong& r1, const slong& c1, const slong& r2, const slong& c2)
{
    window->entries = NULL;

    window->rows = flint_malloc((r2 - r1) * sizeof(acb_ptr));

    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < r2 - r1; i++)
        window->rows[i] = mat->rows[r1 + i] + c1;
        
        flint_cleanup();
    }

    window->r = r2 - r1;
    window->c = c2 - c1;
};


void
arb_mat_add_omp(arb_mat_t& C, const arb_mat_t& A, const arb_mat_t& B, const slong& prec, const slong& num_threads)
{
    if ( (arb_mat_nrows(B) != arb_mat_nrows(A)) ||
         (arb_mat_ncols(B) != arb_mat_ncols(A)) )
    {
        flint_printf("arb_mat_add_omp: incompatible dimensions\n");
        flint_abort();
    }
    
        
    slong r, c;

    r = arb_mat_nrows(C);
    c = arb_mat_ncols(C);

    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < r; i++)
        for (slong j = 0; j < c; j++)
        arb_add(arb_mat_entry(C,i,j),arb_mat_entry(A,i,j),arb_mat_entry(B,i,j),prec);
        
        flint_cleanup();
    }
};


void
arb_mat_cos_entrywise_omp(arb_mat_t& B, const arb_mat_t& A, const slong& prec, const slong& num_threads)
{
    if ( (arb_mat_nrows(A) != arb_mat_nrows(B)) ||
         (arb_mat_ncols(A) != arb_mat_ncols(B)) )
    {
        flint_printf("arb_mat_cos_entrywise_omp: incompatible dimensions\n");
        flint_abort();
    }
    
        
    slong r, c;

    r = arb_mat_nrows(B);
    c = arb_mat_ncols(B);

    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < r; i++)
        for (slong j = 0; j < c; j++)
        acb_cos(arb_mat_entry(B,i,j),arb_mat_entry(A,i,j),prec);
        
        flint_cleanup();
    }
};


void arb_mat_div_entrywise_omp(arb_mat_t& C, const arb_mat_t& A, const arb_mat_t& B, const slong& prec, const slong& num_threads)
{

    if (arb_mat_nrows(A) != arb_mat_nrows(B) ||
        arb_mat_ncols(A) != arb_mat_ncols(B))
    {
        flint_printf("arb_mat_div_entrywise_omp: incompatible dimensions\n");
        flint_abort();
    }
    
    slong r = arb_mat_nrows(A);
    slong c = arb_mat_ncols(A);
    
    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < r; i++)
        {
            for (slong j = 0; j < c; j++)
            {
                arb_div(arb_mat_entry(C, i, j),
                        arb_mat_entry(A, i, j),
                        arb_mat_entry(B, i, j), prec);
            }
        }
        flint_cleanup();
    }
};


void
arb_mat_init_omp(arb_mat_t& mat, const slong& r, const slong& c, const slong& num_threads)
{
    if (r != 0 && c != 0)
    {
        mat->entries = _arb_vec_init(r * c);
        mat->rows = (arb_ptr *) flint_malloc(r * sizeof(arb_ptr));
        omp_set_num_threads(num_threads);
        #pragma omp parallel
        {
            #pragma omp for
            for (slong i = 0; i < r; i++)
            mat->rows[i] = mat->entries + i * c;
            
            flint_cleanup();
        }
    }
    else
        mat->entries = NULL;

    mat->r = r;
    mat->c = c;
};


void
arb_mat_linspace_omp(arb_mat_t& vec, const arb_t& start_arb, const arb_t& end_arb, const slong& prec, const slong& num_threads)
{
    if ( (arb_mat_nrows(vec) != 1 && arb_mat_ncols(vec) != 1) ||
         (arb_mat_nrows(vec) == 1 && arb_mat_ncols(vec) == 1) )
    {
        flint_printf("arb_mat_linspace_omp: incompatible dimensions\n");
        flint_abort();
    }
    
        
    arb_set(arb_mat_entry(vec,0,0), start_arb);
    
    arb_t spacing; arb_init(spacing);
    arb_sub(spacing,end_arb,start_arb,prec);
    arb_div_si(spacing,spacing,np,prec);
    
    omp_set_num_threads(num_threads);
    if (arb_mat_nrows(vec) != 1)
    {
        slong np = arb_mat_nrows(vec);
        arb_mat_set(arb_mat_entry(vec,np-1,0),end_arb);
        #pragma omp parallel
        {
            #pragma omp for
            for(slong i = 1; i < np-1; i++)
            {
                arb_add(arb_mat_entry(vec,i,0),start_arb,prec);
                arb_addmul_si(arb_mat_entry(vec,i,0),spacing,i,prec);
                
            }
            flint_cleanup();
        }
    }
    else
    {
        slong np = arb_mat_ncols(vec);
        arb_mat_set(arb_mat_entry(vec,0,np-1),end_arb);
        #pragma omp parallel
        {
            #pragma omp for
            for(slong i = 1; i < np-1; i++)
            {
                arb_add(arb_mat_entry(vec,0,i),start_arb,prec);
                arb_addmul_si(arb_mat_entry(vec,0,i),spacing,i,prec);
            }
            flint_cleanup();
        }
    }
    arb_clear(spacing);
};


void arb_mat_mul_entrywise_omp(arb_mat_t& C, const arb_mat_t& A, const arb_mat_t& B, const slong& prec, const slong& num_threads)
{
    slong i, j;

    if (arb_mat_nrows(A) != arb_mat_nrows(B) ||
        arb_mat_ncols(A) != arb_mat_ncols(B))
    {
        flint_printf("arb_mat_mul_entrywise_omp: incompatible dimensions\n");
        flint_abort();
    }
    
    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (i = 0; i < arb_mat_nrows(A); i++)
        {
            for (j = 0; j < arb_mat_ncols(A); j++)
            {
                arb_mul(arb_mat_entry(C, i, j),
                        arb_mat_entry(A, i, j),
                        arb_mat_entry(B, i, j), prec);
            }
        }
        flint_cleanup();
    }
};


void
arb_mat_neg_omp(arb_mat_t& dest, const arb_mat_t& src, const slong& num_threads)
{
    if ( (arb_mat_nrows(dest) != arb_mat_nrows(src)) ||
         (arb_mat_ncols(dest) != arb_mat_ncols(src)) )
    {
        flint_printf("arb_mat_neg_omp: incompatible dimensions\n");
        flint_abort();
    }
    
    slong R, C;

    R = arb_mat_nrows(src);
    C = arb_mat_ncols(src);

    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < R; i++)
        for (slong j = 0; j < C; j++)
        arb_neg(arb_mat_entry(dest, i, j), arb_mat_entry(src, i, j));
        
        flint_cleanup();
    }
};


void
arb_mat_ones_omp(arb_mat_t& mat, const slong& num_threads)
{
    slong R, C;

    R = arb_mat_nrows(mat);
    C = arb_mat_ncols(mat);

    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < R; i++)
        for (slong j = 0; j < C; j++)
        arb_one(arb_mat_entry(mat, i, j));
        
        flint_cleanup();
    }
};


void
arb_mat_scalar_addmul_si_omp(arb_mat_t& B, const arb_mat_t& A, const slong& c, const slong& prec, const slong& num_threads)
{
    if ( (arb_mat_nrows(B) != arb_mat_nrows(A)) ||
         (arb_mat_ncols(B) != arb_mat_ncols(A)) )
    {
        flint_printf("arb_mat_scalar_addmul_si_omp: incompatible dimensions\n");
        flint_abort();
    }
        
    slong R, C;

    R = arb_mat_nrows(B);
    C = arb_mat_ncols(B);

    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < R; i++)
        for (slong j = 0; j < C; j++)
        arb_addmul_si(arb_mat_entry(B,i,j),arb_mat_entry(A,i,j),c,prec);
        
        flint_cleanup();
    }
};


void
arb_mat_scalar_div_arb_omp(arb_mat_t& B, const arb_mat_t& A, const arb_t& c, const slong& prec, const slong& num_threads)
{
    if ( (arb_mat_nrows(B) != arb_mat_nrows(A)) ||
         (arb_mat_ncols(B) != arb_mat_ncols(A)) )
    {
        flint_printf("arb_mat_scalar_div_arb_omp: incompatible dimensions\n");
        flint_abort();
    }
    
        
    slong R, C;

    R = arb_mat_nrows(B);
    C = arb_mat_ncols(B);

    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < R; i++)
        for (slong j = 0; j < C; j++)
        arb_div(arb_mat_entry(B,i,j),arb_mat_entry(A,i,j),c,prec);
        
        flint_cleanup();
    }
};


void
arb_mat_scalar_div_si_omp(arb_mat_t& B, const arb_mat_t& A, const slong& c, const slong& prec, const slong& num_threads)
{
    if ( (arb_mat_nrows(B) != arb_mat_nrows(A)) ||
         (arb_mat_ncols(B) != arb_mat_ncols(A)) )
    {
        flint_printf("arb_mat_scalar_div_si_omp: incompatible dimensions\n");
        flint_abort();
    }
    
        
    slong R, C;

    R = arb_mat_nrows(A);
    C = arb_mat_ncols(A);

    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < R; i++)
        for (slong j = 0; j < C; j++)
        arb_div_si(arb_mat_entry(B,i,j),arb_mat_entry(A,i,j),c,prec);
        
        flint_cleanup();
    }
};


void
arb_mat_scalar_mul_arb_omp(arb_mat_t& B, const arb_mat_t& A, const arb_t& c, const slong& prec, const slong& num_threads)
{
    if ( (arb_mat_nrows(B) != arb_mat_nrows(A)) ||
         (arb_mat_ncols(B) != arb_mat_ncols(A)) )
    {
        flint_printf("arb_mat_scalar_mul_arb_omp: incompatible dimensions\n");
        flint_abort();
    }
    
        
    slong R, C;

    R = arb_mat_nrows(A);
    C = arb_mat_ncols(A);

    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < R; i++)
        for (slong j = 0; j < C; j++)
        arb_mul(arb_mat_entry(B,i,j),arb_mat_entry(A,i,j),c,prec);
        
        flint_cleanup();
    }
};


void
arb_mat_scalar_mul_si_omp(arb_mat_t& B, const arb_mat_t& A, const slong& c, const slong& prec, const slong& num_threads)
{
    if ( (arb_mat_nrows(B) != arb_mat_nrows(A)) ||
         (arb_mat_ncols(B) != arb_mat_ncols(A)) )
    {
        flint_printf("arb_mat_scalar_mul_si_omp: incompatible dimensions\n");
        flint_abort();
    }
    
        
    slong R, C;

    R = arb_mat_nrows(A);
    C = arb_mat_ncols(A);

    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < R; i++)
        for (slong j = 0; j < C; j++)
        arb_mul_si(arb_mat_entry(B,i,j),arb_mat_entry(A,i,j),c,prec);
        
        flint_cleanup();
    }
};


void
arb_mat_set_omp(arb_mat_t& dest, const arb_mat_t& src, const slong& num_threads)
{
    slong R, C;

    R = arb_mat_nrows(dest);
    C = arb_mat_ncols(dest);

    omp_set_num_threads(num_threads);
    
    if (dest != src && arb_mat_ncols(src) != 0)
    {
        #pragma omp parallel
        {
            #pragma omp for
            for (slong i = 0; i < R; i++)
            for (slong j = 0; j < C; j++)
            arb_set(arb_mat_entry(dest, i, j), arb_mat_entry(src, i, j));
            
            flint_cleanup();
        }
    }
};


void
arb_mat_sin_cos_entrywise_omp(arb_mat_t& A, arb_mat_t& B, const arb_mat_t& C, const slong& prec, const slong& num_threads);
{
    if ( (arb_mat_nrows(C) != arb_mat_nrows(A)) ||
         (arb_mat_ncols(C) != arb_mat_ncols(A)) ||
         (arb_mat_nrows(C) != arb_mat_nrows(B)) ||
         (arb_mat_ncols(C) != arb_mat_ncols(B)) )
    {
        flint_printf("arb_mat_sin_cos_entrywise_omp: incompatible dimensions\n");
        flint_abort();
    }
    
        
    slong r, c;

    r = arb_mat_nrows(C);
    c = arb_mat_ncols(C);

    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < r; i++)
        for (slong j = 0; j < c; j++)
        arb_sin_cos(arb_mat_entry(A,i,j),arb_mat_entry(B,i,j),arb_mat_entry(C,i,j),prec);
        
        flint_cleanup();
    }
};


void
arb_mat_sin_entrywise_omp(arb_mat_t& B, const arb_mat_t& A, const slong& prec, const slong& num_threads)
{
    if ( (arb_mat_nrows(A) != arb_mat_nrows(B)) ||
         (arb_mat_ncols(A) != arb_mat_ncols(B)) )
    {
        flint_printf("arb_mat_sin_entrywise_omp: incompatible dimensions\n");
        flint_abort();
    }
    
        
    slong r, c;

    r = arb_mat_nrows(B);
    c = arb_mat_ncols(B);

    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < r; i++)
        for (slong j = 0; j < c; j++)
        acb_sin(arb_mat_entry(B,i,j),arb_mat_entry(A,i,j),prec);
        
        flint_cleanup();
    }
};


void
arb_mat_sqrt_entrywise_omp(arb_mat_t& dest, const arb_mat_t& src, const slong& prec, const slong& num_threads)
{
    if ( (arb_mat_nrows(dest) != arb_mat_nrows(src)) ||
         (arb_mat_ncols(dest) != arb_mat_ncols(src)) )
    {
        flint_printf("arb_mat_sqrt_entrywise_omp: incompatible dimensions\n");
        flint_abort();
    }
    
    slong R, C;

    R = arb_mat_nrows(src);
    C = arb_mat_ncols(src);

    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < R; i++)
        for (slong j = 0; j < C; j++)
        arb_sqrt(arb_mat_entry(dest, i, j), arb_mat_entry(src, i, j),prec);
        
        flint_cleanup();
    }
};


void
arb_mat_sub_omp(arb_mat_t& C, const arb_mat_t& A, const arb_mat_t& B, const slong& prec, const slong& num_threads)
{
    if ( (arb_mat_nrows(B) != arb_mat_nrows(A)) ||
         (arb_mat_ncols(B) != arb_mat_ncols(A)) )
    {
        flint_printf("arb_mat_sub_omp: incompatible dimensions\n");
        flint_abort();
    }
    
        
    slong r, c;

    r = arb_mat_nrows(C);
    c = arb_mat_ncols(C);

    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < r; i++)
        for (slong j = 0; j < c; j++)
        arb_sub(arb_mat_entry(C,i,j),arb_mat_entry(A,i,j),arb_mat_entry(B,i,j),prec);
        
        flint_cleanup();
    }
};


void
arb_mat_window_init(arb_mat_t& window, const arb_mat_t& mat, const slong& r1, const slong& c1, const slong& r2, const slong& c2)
{
    window->entries = NULL;

    window->rows = flint_malloc((r2 - r1) * sizeof(arb_ptr));

    #pragma omp parallel
    {
        #pragma omp for
        for (slong i = 0; i < r2 - r1; i++)
        window->rows[i] = mat->rows[r1 + i] + c1;
        
        flint_cleanup();
    }

    window->r = r2 - r1;
    window->c = c2 - c1;
};
