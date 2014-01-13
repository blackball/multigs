#include "model.h"
#include "utils.h"
#include <stdlib.h>

struct model_t {
        int type;
};

struct model_t *
model_new(int type) {

        return 0;
}

void
model_free(struct model_t *model) {
        
}

int
model_fit(struct model_t *model, const struct point_t *ps0, const struct point_t *ps1, int n) {
        
        return 0;
}

int
model_deg(struct model_t *model, const struct point_t *ps0, const struct point_t *ps1, int n) {
        
        return 0;
}

int
model_res(struct model_t *model, const struct point_t *ps0, const struct point_t *ps1, int n, double *out) {
        
        return 0;
}
