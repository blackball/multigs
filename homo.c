#include "utils.h"
#include "svd.h" 
#include "model.h" // struct point_t

static int
homo_deg(const struct point_t *pt0, const struct point_t *pt1, int n) {

        return 0;
}

static int
homo_fit(const struct point_t *pt0, const struct point_t *pt1, int n, struct svdata_t *usv, double model[9]) {
        int i, mini;
        double *u = usv->u, mins;
        const double *s = usv->s;
        const double *v = usv->v;

        if (n != 4) {
                // only support 4 points now.
                retunr -1;
        }
        
#define _ASSIGN(i)                                                      \
        u[i*18+0] = 0;                      u[i*18+1] = 0;                      u[i*18+2] = 0; \
        u[i*18+3] = -pt0[i].x;              u[i*18+4] = -pt0[i].y;              u[i*18+5] = -1; \
        u[i*18+6] = pt1[i].y * pt0[i].x;    u[i*18+7] = pt1[i].y * pt0[i].y;    u[i*18+8] = pt1[i].y; \
        u[i*18+9] = pt0[i].x;               u[i*18+10] = pt0[i].y;              u[i*18+11] = 1; \
        u[i*18+12] = 0;                     u[i*18+13] = 0;                     u[i*18+14] = 0; \
        u[i*18+15] = -pt1[i].x * pt0[i].x;  u[i*18+16] = -pt1[i].x * pt0[i].y;  u[i*18+17] = -pt1[i].x
 
        _ASSIGN(0);
        _ASSIGN(1);
        _ASSIGN(2);
        _ASSIGN(3);
        
#undef _ASSIGN
 
        if (!svd(usv)) {
                return -1;
        }
 
        mini = 0;
        mins = s[0];
        for (i = 1; i < 9; ++i) {
                if (s[i] < mins) {
                        mins = s[i];
                        mini = i;
                }
        }

        for (i = 0; i < 9; ++i) {
                model[i] = v[mini * 9 + i];
        }
        
        return 0;
}

static int
homo_res(const double H[9], const struct point_t *pts0, const struct point_t *pts1, int n, double *out) {
        int i = 0;
        for (; i < n ; ++i) {
                out[i] = .0;
        }
        return 0;
}
