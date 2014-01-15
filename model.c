#include "model.h"
#include "utils.h"
#include "svd.h"
#include <stdlib.h>

#include "norm.c" /* implementation of points normalization for DLT */
#include "homo.c" /* implementation of 4-point homography estimation */

struct model_t {
        int type;
        double H[9]; // model
        struct svdata_t *svdata89; /* size of u is 8x9 */
};

struct model_t *
model_new(int type) {
        struct model_t *model = malloc(sizeof(*model));
        DASSERT(model);
        model->type = type;
        if (type == MODEL_TYPE_HOMO) {
                model->svdata89 = svdata_new(8, 9);
        }
        /* add more type here! */
        else {
                model_free(model);
                model = 0;
                DASSERT(0); // not supported;
        }
        
        return model;
}

void
model_free(struct model_t *model) {
        if (model) {
                svdata_free(model->svdata89);
                free(model);
        }
}

int
model_fit(struct model_t *model, const struct point_t *ps0, const struct point_t *ps1, int n) {
        if (model->type == MODEL_TYPE_HOMO) {
                struct point_t nps0[4], nps1[4];
                double T0[3], T1[3];
                
                DASSERT(n == 4);
#if MODEL_NEED_NORMALIZE == 1
                norm_get_condition(ps0, n, T0);
                norm_get_condition(ps1, n, T1);
                norm_transform(x, ps0, nps0, 4);
                norm_transform(x, ps1, nps1, 4);
#else
                nps0[0] = ps0[0];
                nps0[1] = ps0[1];
                nps0[2] = ps0[2];
                nps0[3] = ps0[3];
                nps1[0] = ps1[0];
                nps1[1] = ps1[1];
                nps1[2] = ps1[2];
                nps1[3] = ps1[3];
#endif
                return homo_fit(nps0, nps1, n, model->svdata89, model->H);
        }
        /* add others here */
        else {
                DASSERT(0);
                return 1;
        }
        return 0;
}

int
model_deg(struct model_t *model, const struct point_t *ps0, const struct point_t *ps1, int n) {
        if (model->type == MODEL_TYPE_HOMO) {
                return homo_deg(ps0, ps1, n);
        }
        /* add others here */
        else {
                DASSERT(0);
                return 1;
        }
        return 0;
}

int
model_res(struct model_t *model, const struct point_t *ps0, const struct point_t *ps1, int n, double *out) {
        if (model->type == MODEL_TYPE_HOMO) {
                
        }
        /* add others here */
        else {
                DASSERT(0);
                return -1;
        }
        return 0;
}
