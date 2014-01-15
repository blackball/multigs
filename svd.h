#ifndef MULTIGS_SVD_H
#define MULTIGS_SVD_H


struct svdata_t {
        int m, n;
        double *u;  /* mxn */
        double *s;  /* 1xn */
        double *v;  /* nxn */
        double *rv; /* 1xn, buffer */       
};

/* Create u,s,v that used in SVD, *u* is a mxn martix */
struct svdata_t svdata_new(int m, int n);

void svdata_free(struct svdata_t *svdata);

/* SVD, return 0 if OK, else -1 */
int svd(struct svdata_t *svdata);

#endif
