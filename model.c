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
