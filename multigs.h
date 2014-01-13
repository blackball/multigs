#ifndef MULTIGS_H
#define MULTIGS_H

struct point_t; // just used in model estimation algorithms
struct multigs_t;
struct multigs_t * multigs_new(int hypo_size, int init_size);
void multigs_free(struct multigs_t *gs);
int multigs_sampling(struct multigs_t *gs, const struct point_t *ps0, const struct point_t *ps1);

#endif
