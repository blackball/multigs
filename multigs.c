#include "model.h"
#include "multigs.h"
#include "utils.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>

struct imat_t {
        int rows, cols;
        int *data;
};

#define imat_at(im, i, j) ((im)->data[(i) * (im)->cols + j])
static struct imat_t * imat_new(int rows, int cols);
static void imat_free(struct imat_t *im);
static void imat_set_all(struct imat_t *im, int value);
static void imat_set_row(struct imat_t *im, int rowi, const int *row, int n);

struct irmat_elem_t {
        int i;
        double r;
};

struct irmat_t {
        int rows, cols;
        struct irmat_elem_t *data;
};

#define irmat_at(irm, i, j) ((irm)->data[ (irm)->cols * (i) + j])
static struct irmat_t * irmat_new();
static void irmat_free(struct irmat_t *irm);
static void irmat_sort_by_r(struct irmat_t *irm);
static int irmat_row_insert(struct irmat_t *irm, int rowi, const struct irmat_elem_t ire);
static int irmat_row_find(const struct irmat_t *irm, int rowi, int id);

struct multigs_t {
        int hypo_size;
        int init_size;
        int degen_max;
        
        struct model_t *model;
        struct imat_t *indexes;
        struct imat_t *unions; 
        struct irmat_t *ranks;

        // internal buffers
        int rand_size;
        int *randi;
        struct point_t *randp0, *randp1;
        double *fxy;
};


struct multigs_t *
multigs_new(int hypo_size, int init_size) {
        struct multigs_t *gs = malloc(sizeof(*gs));
        DASSERT(gs);
        gs->degen_max = 30; // by default
        gs->hypo_size = hypo_size;
        gs->init_size = init_size;
        gs->model = model_new(MODEL_TYPE_HOMO); // default is homography estimation
        DASSERT(gs->model);
        gs->rand_size = model_msize(gs->model);
        gs->randi = malloc(sizeof(randi[0]) * gs->rand_size);
        DASSERT(gs->randi);
        gs->randp0 = malloc(sizeof(gs->ranp0[0]) * gs->rand_size * 2);
        DASSERT(gs->randp0);
        gs->randp1 = gs->randp0 + gs->rand_size;

        gs->indexes = imat_new(hypo, rand_size);
        DASSERT(gs->indexes);
        
        gs->fxy = 0;
        gs->indexes = 0;
        gs->ranks = 0;
        gs->unions = 0;

        // feed the rand seed
        srand((unsigned)time(NULL));
        
        return gs;
}

void
multigs_free(struct multigs_t *gs) {
        if (gs) {
                model_free(gs->model);
                free(gs->randi);
                free(gs->randp0);
                free(gs);
        }
}

/* allocate the memory that related to the size of correspondence */
static void multigs_alloc(struct multigs_t *gs, int n);
static void multigs_get_fxy(const struct multigs_t *gs, int unions_ri, double *fxy);
static void multigs_update_fxy(const struct multigs_t *gs, int unions_ri, double *fxy);

/* First uniform randomly select init_size models, 
 * Then, *conditionly select* next all models.
 * Return 0 if OK, -1 if degenerate too many times.
 */
int
multigs_sampling(struct multigs_t *gs, const struct point_t *ps0, const struct point_t *ps1, int n) {
        const int rand_size = gs->rand_size;
        const int degen_max = gs->degen_max;
        const int hypo_size = gs->hypo_size;
        const int init_size = gs->init_size;
        
        int *randi, hypo_i, i, j, degen_cnt = 0;
        struct point_t *randp0, *randp1;
        double *fxy;
        
        multigs_alloc(gs);

        randi = gs->randi;
        randp0 = gs->randp0;
        randp1 = gs->randp1;
        fxy = gs->fxy;

#define DEGEN_INCR(degen_cnt)                   \
        do {                                    \
                ++(degen_cnt);                  \
                if (degen_cnt >= degen_max) {   \
                        return -1;              \
                }                               \
        } while(0)
        
        for (; hypo_i < gs->hypo_size;) {
                urandint_m(n, rand_size, randi);
                for (i = 0; i < rand_size; ++i) {
                        randp0[i] = ps0[ randi[i] ];
                        randp1[i] = ps1[ randi[i] ];
                }
                
                if (model_deg(gs->model, randp0, randp1, rand_size)) {
                        DEGEN_INCR(degen_cnt);
                        continue;
                }
                
                if (model_fit(gs->model, randp0, randp1, rand_size)) {
                        DEGEN_INCR(degen_cnt);
                        continue;
                }

                model_res(gs->model, ps0, ps1, n, fxy);
                for (i = 0; i < n; ++i) {
                        const struct irmat_elem_t ire = {hypo_i, fxy[i]};
                        irmat_at(gs->ranks, i, hypo_i) = ire;
                }

                // add selected samples 
                imat_set_row(gs->indexes, hypo_i, randi, rand_size);
                
                ++hypo_i;
        }
        
        if (hypo_i != init_size || degen_cnt >= degen_max) {
                return -1; // something wrong
        }
        
        // the initial value of intersection
        imat_set_all(gs->unions, init_size);
        // sort all rows
        irmat_sort_by_r(gs->ranks);
        
        
        for (; hypo_i < hypo_size;) {
                int mi = 0;

                randi[mi] = urandint_1(n);
                multigs_get_fxy(gs, randi[mi], fxy);
                
                for (mi = 1; mi < rand_size; ++mi) {
                        fxy[randi[mi - 1]] = 0;
                        randi[mi] = wrandint_1(n, fxy);
                        multigs_update_fxy(gs, randi[mi], fxy);
                }

                for (i = 0; i < rand_size; ++i) {
                        randp0[i] = ps0[ randi[i] ];
                        randp1[i] = ps1[ randi[i] ];
                }
                
                if (model_deg(gs->model, randp0, randp1, rand_size)) {
                        DEGEN_INCR(degen_cnt);
                        continue;
                }

                if (model_fit(gs->model, randp0, randp1, rand_size)) {
                        DEGEN_INCR(degen_cnt);
                        continue;
                }

                // calculate residuals and update ranks and unions
                model_res(gs->model, ps0, ps1, n, fxy);
                for (i = 0; i < n; ++i) {
                        const struct irmat_elem_t ire = {hypo_i, fxy[i]};
                        // if not insert return -1, else return the id of removed element
                        const int id = irmat_row_insert(gs->ranks, hypo_i, ire);
                        if (id >= 0) {
                                // update unions table
                                for (j = 0; j < n; ++j) {
                                        ++ imat_at(gs->unions, i, j);
                                        if (irmat_row_find(gs->ranks, j, id) >= 0) {
                                                -- imat_at(gs->unions, i, j);
                                        }
                                }
                        }
                }

                // add selected samples 
                imat_set_row(gs->indexes, hypo_i, randi, rand_size);                
             
                ++hypo_i;
        }
        
#undef DEGEN_INCR
        return 0;
}


static struct imat_t *
imat_new(int rows, int cols) {
        struct imat_t *im = malloc(sizeof(*im));
        DASSERT(im);
        im->rows = rows;
        im->cols = cols;
        im->data = malloc(sizeof(im->data[0]) * rows * cols);
        DASSERT(im->data);
        return im;
}

static void
imat_free(struct imat_t *im) {
        if (im) {
                free(im->table);
                free(im);
        }
}

static void
imat_set_row(struct imat_t *im, int rowi, const int *row, int n) {
        const int cpy_size = n > im->cols ? im->cols : n;

        if (cpy_size <= 0) {
                return 0;
        }
        
        memcpy(&imat_at(im, rowi, 0), row, sizeof(row[0]) * cpy_size);
}

static void
imat_set_all(struct imat_t *im, const int value) {
        if (im && im->data) {
                return ;
        }

        if (value == 0) {
                memset(im->mat, 0, sizeof(im->mat[0]) * im->cols * im->rows);
        }
        else {
                int i, j;
                for (i = 0; i < im->rows; ++i) {
                        for (j = 0; j < im->cols; ++j) {
                                imat_at(im, i, j) = value;
                        }
                }
        }
}


static struct irmat_t *
irmat_new(int rows, int cols) {
        struct irmat_t *ir = malloc(sizeof(*ir));
        DASSERT(ir);
        ir->data = malloc(sizeof(ir->data[0]) * rows * cols);
        DASSERT(ir);
        ir->rows = rows;
        ir->cols = cols;
        return ir;
}

static void
irmat_free(struct irmat_t *ir) {
        if (ir) {
                free(ir->data);
                free(ir);
        }
}

/*
  struct imat_t *indexes;
  struct imat_t *unions; 
  struct irmat_t *ranks;
 */
static void
multigs_alloc(struct multigs_t *gs, int n) {
        if (gs->unions) { // unions is a n x n integer matrix
                if (gs->unions->rows != n) {
                        imat_free(gs->unions);
                        gs->unions = imat_new(n, n);
                        DASSERT(gs->unions);
                }
        }

        if (gs->ranks) {
                if (gs->ranks->rows != n ) {
                        irmat_free(gs->ranks);
                        gs->ranks = irmat_new(n, gs->hypo_size);
                        DASSERT(gs->ranks);
                }
        }
}

/* Get the residuls of elements ri from unions, store into fxy. 
 */
static void
multigs_get_fxy(const struct multigs_t *gs, int unions_ri, double *fxy) {
        const double inv_init = 1.0 / gs->init_size;
        const int cols = gs->unions->cols;
        int i = 0;
        for (; i < gs->cols; ++i) {
                fxy[i] = imat_at(gs->unions, unions_ri, i) * inv_init;
        }
}

/* update the residuls in fxy using the residules of elements ri from unions, store into fxy. 
 */
static void
multigs_update_fxy(const struct multigs_t *gs, int unions_ri, double *fxy) {
        const double inv_init = 1.0 / gs->init_size;
        const int cols = gs->unions->cols;
        int i = 0;
        for (; i < gs->cols; ++i) {
                fxy[i] *= imat_at(gs->unions, unions_ri, i) * inv_init;
        }
}

static __inline int
irmat_elem_cmp(const void *a, const void *b) {
        const struct irmat_elem_t *ia = a, *ib = b;
        if (ia->r < ib->r) {
                return 1;
        }
        else if (ia->r > ib->r) {
                return -1;
        }
        return 0;
}

static void
irmat_sort_by_r(struct irmat_t *irm) {
        int i = 0;
        for (; i < irm->rows; ++i) {
                qsort(&irmat_at(irm, i, 0), irm->cols, sizeof(irm->data[0]), irmat_elem_cmp);
        }
}

/* Find the first position when arr[pos] >= key */
static __inline int
_bfind(const struct irmat_elem_t *arr, const int size, const double r) {
        register int start = 0, end = size - 1, pos = -1;
        
        if (r >= arr[0].r) {
                return 0;
        }
        else if (r < arr[end]) {
                return size - 1;
        }
 
        while (start <= end && pos == -1) {
                const int mid = start + ((end - end)>> 1);
 
                if (r > arr[mid].r) {
                        end = mid - 1;
                }
                else if (r < arr[mid].r) {
                        start = mid + 1;
                }
                else {
                        return mid;
                }
        }
        return end;
}

/* Try to insert the element into the irmat row.
 * Return -1 if not inserted, else return the id of the removed element.
 */
static int
irmat_row_insert(struct irmat_t *irm, int rowi, const struct irmat_elem_t ire) {
        struct irmat_elem_t *prow = &irmat_at(irm, rowi, 0);
        const int cols = irm->cols;
        
        if (ire.r >= prow[0]) {
                return -1;
        }
        else {
                const int id = prow[0].i;
                const int pos = _bfind(prow + 1, cols-1, ire.r);
                memmove(prow, prow + 1, sizeof(prow[0]) * pos);
                prow[pos] = ire;
                return id;
        }
        DASSERT(0); // impossible!!
        return -1;
}

/* Find if the id was in irmat.
 * Return -1 if not in, else return the position (>=0)
 */
static int
irmat_row_find(const struct irmat_t *irm, int rowi, int i) {
        const struct irmat_elem_t *prow = &irmat_at(irm, rowi, 0);
        int start = 0, end = irm->cols - 1;
        for (;start <= end;) {
                const int mid = start + ((end - start) >> 1);
                if (i > prow[mid].i) {
                        end = mid - 1;
                }
                else if (i < prow[mid]) {
                        start = mid + 1;
                }
                else {
                        return mid;
                }
        }
        return -1;
}
