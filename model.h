#ifndef MULTIGS_MODEL_H
#define MULTIGS_MODEL_H

#define MODEL_TYPE_HOMO 0
#define MODEL_TYPE_FUND 1

struct point_t {
        double x[2];
};
        
struct model_t;

struct model_t * model_new(int type);
void model_free(struct model_t *model);

/* Model estimation. 
 * Return 0 if success, else -1
 */
int model_fit(struct model_t *model, const struct point_t *ps0, const struct point_t *ps1, int n);

/* Check if the selected correspondences are not good.
 * Return 1 if degenerate, else 0;
 */
int model_deg(struct model_t *model, const struct point_t *ps0, const struct point_t *ps1, int n);

/* Calculate residuals.
 * The matrix in model should be calculated before, else the results will be random...
 * Return 0, at the time.
 */
int model_res(struct model_t *model, const struct point_t *ps0, const struct point_t *ps1, int n, double *out);

#endif
