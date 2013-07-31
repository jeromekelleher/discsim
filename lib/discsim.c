/*
** Copyright (C) 2013 Jerome Kelleher <jerome.kelleher@ed.ac.uk>
**  
** This file is part of ercs.
** 
** ercs is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
** 
** ercs is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with ercs.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <sys/times.h>
#include <unistd.h>
#include <stdint.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_histogram.h>

static void
fatal_error(const char *message)
{
    fprintf(stderr, "%s\n", message);
    exit(1);
}

void *
xmalloc(size_t size)
{
    register void *value = malloc(size);
    if (value == NULL) {
        fatal_error("virtual memory exhausted");
    }
    return value;
}

void *
xcalloc(size_t count, size_t eltsize)
{
    register void *value = calloc(count, eltsize);
    if (value == NULL) {
        fatal_error("virtual memory exhausted");
        abort();
    }
    return value;
}

static int 
avl_set_compare(const void *pa, const void *pb)
{
    uint64_t a = *((uint64_t *) pa);
    uint64_t b = *((uint64_t *) pb);
    return (a > b) - (a < b);
}

static int 
avl_set_map_compare(const void *pa, const void *pb)
{
    int a = ((set_map_value_t *) pa)->key;
    int b = ((set_map_value_t *) pb)->key;
    return (a > b) - (a < b);
}

static int 
avl_int_map_compare(const void *pa, const void *pb)
{
    int a = ((int_map_value_t *) pa)->key;
    int b = ((int_map_value_t *) pb)->key;
    return (a > b) - (a < b);
}


/* 
 * Returns the expected density of ancestors at equilibrium in a simulation 
 * with nu parents and impact u.
 */
double 
equilibrium_density(double nu, double u)
{
    double N = nu / u;
    double lu = log(1.0 - u);
    double numer = gsl_sf_lambert_W0(N * pow(1 - u, N) * lu);
    double ret = (N  - numer / lu) / M_PI;
    return ret;
}


/*
 * Returns the square of the distance between the two specified points on a 
 * square torus of side R.
 */
static inline double
torus_squared_distance(const double *p1, const double *p2, double L) 
{
    double xabs = fabs(p2[0] - p1[0]);
    double yabs = fabs(p2[1] - p1[1]);
    double xd = GSL_MIN(xabs, L - xabs);
    double yd = GSL_MIN(yabs, L - yabs);
    return xd * xd + yd * yd;
}

/*
 * Returns the area of the intersection of two discs of radius r with centres
 * distance d apart.
 */
static inline double 
disc_intersection_area(double d, double r)
{
    double ret = 0.0;
    if (d < 2 * r) {
        ret = 2 * r * r * acos(d / (2 * r)) - d * sqrt(4 * r * r - d * d) / 2;
    }
    return ret;
}

/*
 * Uses the specified random variate u to return the index of the specified 
 * list with the appropriate probability.
 */
int 
probability_list_select(double *probabilities, const unsigned int n, 
        const double u) 
{
    int ret = 0;
    int x;
    unsigned int i;
    double lhs, rhs;
    if (n == 0) {
        ret = -1;
    } else if (n == 1) {
        ret = 0;
    } else {
        x = -1;
        i = 0;
        lhs = 0.0; 
        rhs = probabilities[0];
        while (x == -1) {
            if (lhs <= u && u < rhs) {
                x = (int) i;
            } else {
                i++;
                lhs = rhs;
                rhs += probabilities[i];
            }
        }
        ret = x; 
    }
    return ret;
}

/*
 * Updates the specified point such that it contains a point uniformly distributed
 * within a circle of radius r, centered at the point centre on a square torus of 
 * side L.
 *
 * This uses the canonical method for generating a point uniformly distributed within
 * the unit circle (from TAOCP), then scales and takes the modulus appropriately.
 */
static inline void 
random_point_torus_disc(double *p, const double *centre, const double r, 
        const double L, gsl_rng *generator)
{
    register double s = 1.1;
    register double x;
    register double y;
    while (s >= 1) {   
        x = 2 * gsl_rng_uniform(generator) - 1;
        y = 2 * gsl_rng_uniform(generator) - 1;
        s = (x * x) + (y * y);
    }
    x = centre[0] + (x * r);
    y = centre[1] + (y * r);
    p[0] = fmod(x + L, L);
    p[1] = fmod(y + L, L); 
}

static inline int 
pixel_coord_to_index(unsigned int x, unsigned int y, unsigned int N)  
{
    int ret = x * N + y;
    return ret;
}

static inline void 
index_to_pixel_coord(int index, unsigned int N, unsigned int *v)
{
    v[0] = index / N;
    v[1] = index % N;
}

double
beta(int n, int k, double u)
{
    double ret = gsl_sf_choose(n, k);
    ret *= gsl_pow_int(u, k) * gsl_pow_int(1.0 - u, n - k);
    ret /= 1.0 - gsl_pow_int(1.0 - u, n);
    //printf("beta %d %d %f = %.14f\n", n, k, u, ret);
    return ret;

}


void
sim_print_parameters(sim_t *self)
{
    unsigned int j;
    double *x;
    printf("# L = %f\n", self->L);
    printf("# m = %u\n", self->m);
    printf("# s = %f\n", self->s);
    printf("# u = %f\n", self->u);
    printf("# r = %f\n", self->r);
    printf("# nu = %u\n", self->nu);
    printf("# rho = %f\n", self->rho);
    printf("# lambda = %f\n", self->lambda);
    printf("# X = ");
    for (j = 0; j < self->n; j++) {
        x = self->X + 2 * j;
        printf("(%f, %f) ", x[0], x[1]); 
    }
    printf("\n# N = %f\n", self->nu / self->u);
}

static inline double 
clamp(double x, double min, double max)
{
    double ret = x;
    if (x < min) {
        ret = min;
    } else if (x > max) {
        ret = max;
    }
    return ret;
}


/* Returns 1 if a disc of radius r intersects with a pixel of side s
 * with bottom-left corner x*s, y*s.
 */
static int 
disc_intersects_pixel(double *z, int x, int y, double r, double s)
{
    // Get the closest point in the pixel to the centre of the disc
    double closest_x = clamp(z[0], x * s, (x + 1) * s);
    double closest_y = clamp(z[1], y * s, (y + 1) * s);
    // Calculate the distance between the disc's center and this closest point
    double dx2 = gsl_pow_2(z[0] - closest_x);
    double dy2 = gsl_pow_2(z[1] - closest_y);
    return  dx2 + dy2 < gsl_pow_2(r);
}


/* Gets the pixels for the disc at the specified location and updates the pixel 
 * buffer.
 */
static void 
get_pixels_general(double r, double s, int N, double *z, unsigned int *p)
{
    int x_min = (int) floor((z[0] - r) / s);
    int x_max = (int) ceil((z[0] + r) / s);
    int y_min = (int) floor((z[1] - r) / s);
    int y_max = (int) ceil((z[1] + r) / s);
    int x, y;
    int coord[2];
    p[0] = 0;
    for (y = y_min; y <= y_max; y++) {
        for (x = x_min; x <= x_max; x++) {
            if (disc_intersects_pixel(z, x, y, r, s)) {
                coord[0] = (N + x) % N;
                coord[1] = (N + y) % N;
                //printf("generated: %d %d\n", coord[0], coord[1]);
                p[0]++; 
                p[p[0]] = pixel_coord_to_index(coord[0], coord[1], N);
            }
        }
    } 
}

static void
get_pixels(double r, double s, int N, double *z, unsigned int *p)
{
    double r2 = r * r;
    double a[2], b[2], c[2], d[2], x[2];
    double L = s * N;
    int v[2];
    assert(s >= 1.0);
    v[0] = (unsigned int) floor(z[0] / s); 
    v[1] = (unsigned int) floor(z[1] / s); 
    p[0] = 1;
    p[p[0]] = pixel_coord_to_index(v[0], v[1], N);
    a[0] = s * v[0];
    a[1] = s * v[1];
    b[0] = a[0] + r;
    b[1] = a[1] + r;
    d[0] = a[0] + s;
    d[1] = a[1] + s;
    c[0] = d[0] - r;
    c[1] = d[1] - r;
    if (z[0] < b[0]) {
        p[0]++;
        p[p[0]] = pixel_coord_to_index((v[0] - 1 + N) % N, v[1], N);
    } 
    if (z[0] > c[0]) {
        p[0]++;
        p[p[0]] = pixel_coord_to_index((v[0] + 1) % N, v[1], N);
    }
    if (z[1] < b[1]) {
        p[0]++;
        p[p[0]] = pixel_coord_to_index(v[0], (v[1] - 1 + N) % N, N);
    }
    if (z[1] > c[1]) {
        p[0]++;
        p[p[0]] = pixel_coord_to_index(v[0], (v[1] + 1) % N, N);
    }
    /* check the corners clocwise from v */ 
    if (z[0] < b[0] && z[1] < b[1]) {
        if (torus_squared_distance(a, z, L) < r2) {
            p[0]++;
            p[p[0]] = pixel_coord_to_index((v[0] - 1 + N ) % N, 
                    (v[1] - 1 + N) % N, N);
        } 
    } 
    if (z[0] < b[0] && z[1] > c[1]) {
        x[0] = a[0];
        x[1] = d[1];
        if (torus_squared_distance(x, z, L) < r2) {
            p[0]++;
            p[p[0]] = pixel_coord_to_index((v[0] - 1 + N) % N, 
                    (v[1] + 1) % N, N);
        }
    } 
    if (z[0] > c[0] && z[1] > c[1]) {
        if (torus_squared_distance(d, z, L) < r2) {
            p[0]++;
            p[p[0]] = pixel_coord_to_index((v[0] + 1) % N, 
                    (v[1] + 1) % N, N);
        }
    } 
    if (z[0] > c[0] && z[1] < b[1]) {
        x[0] = d[0];
        x[1] = a[1];
        if (torus_squared_distance(x, z, L) < r2) {
            p[0]++;
            p[p[0]] = pixel_coord_to_index((v[0] + 1) % N, 
                    (v[1] - 1 + N) % N, N);
        }
    }
#ifdef CHECK_PIXELS 
    unsigned int j;
    unsigned int p2[10];
    get_pixels_general(r, s, N, z, p2);
    assert(p[0] == p2[0]);
    gsl_sort_uint(p + 1, 1, p[0]);
    gsl_sort_uint(p2 + 1, 1, p2[0]);
    for (j = 1; j <= p2[0]; j++) {
        assert(p[j] == p2[j]);
    }
#endif

}


/* Gets the pixels for the disc at the specified location and updates the pixel 
 * buffer.
 */
static void 
sim_get_disc_pixels(sim_t *self, double *z)
{
    if (self->s < 1) {
        get_pixels_general(self->r, self->s, self->N, z, self->pixel_buffer);
    } else {
        get_pixels(self->r, self->s, self->N, z, self->pixel_buffer);
    }

}


/*
 * Sets the max occupancy to reflect the parameters of the simulation 
 * using the equilibrium prediction, using the specified fraction of headroom
 * over this estimate.
 */
static void
sim_set_max_occupancy(sim_t *self, double headroom)
{
    double delta = equilibrium_density(self->nu, self->u);
    double o = delta * (gsl_pow_2(self->s) + 4 * self->s + M_PI);
    o += o * headroom;
    self->max_occupancy = (unsigned int) o; 

}

static void
sim_alloc(sim_t *self)
{
    unsigned int j, k, max;
    double *b;
    self->N = self->L / self->s;
    if (fmod(self->L, self->s) != 0.0) {
        fatal_error("L must be a multiple of s");
    }
    self->rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(self->rng, self->random_seed);
    self->P = xmalloc(gsl_pow_2(self->N) * sizeof(avl_tree_t));
    /* set up buffers */
    self->max_disc_pixels = gsl_pow_2(4 * ((int)(self->r / self->s) + 1));
    self->pixel_buffer = xmalloc(self->max_disc_pixels * sizeof(unsigned int));
    self->probability_buffer = xmalloc(self->max_occupancy * sizeof(double));
    self->intersected_buffer = xmalloc(
            (self->max_occupancy + 1) * sizeof(unsigned int));
    self->child_buffer = xmalloc(
            (self->max_occupancy + 1) * sizeof(unsigned int));
    self->parent_buffer = xmalloc((self->nu + 1) * sizeof(avl_node_t *));
    /* Set up the AVL trees */
    avl_init_tree(&self->Q, avl_set_map_compare, NULL);
    for (j = 0; j < gsl_pow_2(self->N); j++) {
        avl_init_tree(&self->P[j], avl_set_compare, NULL);
        
    }   
    /* Set up the set memory */
    max = self->max_disc_pixels * self->max_sample_size;
    self->avl_set_value_mem = xmalloc(max * sizeof(uint64_t));
    self->avl_set_node_mem = xmalloc(max * sizeof(avl_node_t));
    self->avl_set_node_heap = xmalloc(max * sizeof(avl_node_t *));
    for (j = 0; j < max; j++) {
        self->avl_set_node_heap[j] = self->avl_set_node_mem + j;
        avl_init_node(self->avl_set_node_heap[j], 
                &self->avl_set_value_mem[j]);
    }
    self->avl_set_node_heap_top = max - 1; 
    /* Set up the point memory */
    max = self->max_sample_size;
    self->individual_mem = xmalloc(max * sizeof(individual_t));
    self->individual_heap = xmalloc(max * sizeof(individual_t *));
    for (j = 0; j < max; j++) {
        self->individual_heap[j] = self->individual_mem + j;
    }
    self->individual_heap_top = max - 1; 
    /* Set up the set map memory */
    max = self->max_occupancy;
    self->avl_set_map_value_mem = xmalloc(max * sizeof(set_map_value_t));
    self->avl_set_map_node_mem = xmalloc(max * sizeof(avl_node_t));
    self->avl_set_map_node_heap = xmalloc(max * sizeof(avl_node_t *));
    for (j = 0; j < max; j++) {
        self->avl_set_map_node_heap[j] = self->avl_set_map_node_mem + j;
        avl_init_node(self->avl_set_map_node_heap[j], 
                &self->avl_set_map_value_mem[j]);
    }
    self->avl_set_map_node_heap_top = max - 1; 
    /* Set up the int map memory */
    max = (self->n * gsl_max(self->nu, 2)) * self->m;
    self->avl_int_map_value_mem = xmalloc(max * sizeof(int_map_value_t));
    self->avl_int_map_node_mem = xmalloc(max * sizeof(avl_node_t));
    self->avl_int_map_node_heap = xmalloc(max * sizeof(avl_node_t *));
    for (j = 0; j < max; j++) {
        self->avl_int_map_node_heap[j] = self->avl_int_map_node_mem + j;
        avl_init_node(self->avl_int_map_node_heap[j], 
                &self->avl_int_map_value_mem[j]);
    }
    self->avl_int_map_node_heap_top = max - 1; 
    /* Set up memory for ancestry */
    self->eta = xmalloc(self->m * sizeof(unsigned int));
    self->pi = xmalloc(self->m * sizeof(int *));
    self->tau = xmalloc(self->m * sizeof(double *));
    self->pi_mem = xmalloc(2 * self->n * self->m * sizeof(int));
    self->tau_mem = xmalloc(2 * self->n * self->m * sizeof(double));
    for (j = 0; j < self->m; j++) {
        self->pi[j] = self->pi_mem + j * 2 * self->n; 
        self->tau[j] = self->tau_mem + j * 2 * self->n; 
    }
    self->coalescence_map = xmalloc(self->m * sizeof(unsigned int));
    self->coalesced_loci = xmalloc((self->m + 1) * sizeof(unsigned int));
    /* Precalculate ubar */
    self->ubar = xmalloc((self->max_occupancy + 1) * sizeof(double));
    for (j = 0; j <= self->max_occupancy; j++) {
        self->ubar[j] = 1.0 - gsl_pow_int(1.0 - self->u, j);
    }
    /* Precalculate beta */
    b = xmalloc((self->max_occupancy + 1) * sizeof(double));
    b[0] = 0.0;
    self->beta_distributions = xmalloc((self->max_occupancy + 1) 
            * sizeof(gsl_ran_discrete_t *));
    for (k = 1; k <= self->max_occupancy; k++) {
        /* first figure out how many places we need to get below threshold */
        j = 1;
        while (beta(k, j, self->u) > self->beta_threshold && j < k) {
            j++;
        }
        max = j;
        for (j = 1; j <= max; j++) {
            b[j] = beta(k, j, self->u);
        }
        self->beta_distributions[k] = gsl_ran_discrete_preproc(max + 1, b);
        if (self->beta_distributions[k] == NULL) {
            fatal_error("Out of memory\n");
        }
    }
    free(b);
    /* allocate counter matrices */
    self->total_C_size = xmalloc(gsl_pow_2(self->N) * sizeof(unsigned long));
    self->total_S_size = xmalloc(gsl_pow_2(self->N) * sizeof(unsigned long));
    self->total_G_size = xmalloc(gsl_pow_2(self->N) * sizeof(unsigned long));
    self->total_H_size = xmalloc(gsl_pow_2(self->N) * sizeof(unsigned long));
    self->total_E_size = xmalloc(gsl_pow_2(self->N) * sizeof(unsigned long));
    self->hit_rate_distribution = gsl_histogram_alloc(100);
    gsl_histogram_set_ranges_uniform(self->hit_rate_distribution, 0.0, 100.0);

}

static void
sim_free(sim_t *self)
{
    unsigned int j;
    gsl_rng_free(self->rng);
    free(self->P);
    free(self->pixel_buffer);
    free(self->probability_buffer);
    free(self->intersected_buffer);
    free(self->child_buffer);
    free(self->parent_buffer);
    free(self->avl_set_value_mem);
    free(self->avl_set_node_mem);
    free(self->avl_set_node_heap);
    free(self->individual_mem);
    free(self->individual_heap);
    free(self->avl_set_map_value_mem);
    free(self->avl_set_map_node_mem);
    free(self->avl_set_map_node_heap);
    free(self->avl_int_map_value_mem);
    free(self->avl_int_map_node_mem);
    free(self->avl_int_map_node_heap);
    free(self->pi);
    free(self->tau);
    free(self->eta);
    free(self->pi_mem);
    free(self->tau_mem);
    free(self->coalescence_map);
    free(self->coalesced_loci);
    free(self->ubar);
    for (j = 1; j < self->max_occupancy + 1; j++) {
        gsl_ran_discrete_free(self->beta_distributions[j]);
    }
    free(self->beta_distributions);
    /* counters */
    free(self->total_S_size);
    free(self->total_C_size);
    free(self->total_G_size);
    free(self->total_H_size);
    free(self->total_E_size);
    gsl_histogram_free(self->hit_rate_distribution);
}

/* TODO These function names are inconsistent -  the free functions keep 
 * the avl prefix and the alloc functions drop it.
 */

static avl_node_t *
sim_alloc_set_node(sim_t *self, uint64_t value) 
{
    avl_node_t *node = NULL;
    if (self->avl_set_node_heap_top < 0) {
        fatal_error("Out of set space");
    }
    node = self->avl_set_node_heap[self->avl_set_node_heap_top];
    self->avl_set_node_heap_top--;
    *((uint64_t *) node->item) = value; 
    return node;
}

static void
sim_free_avl_set_node(sim_t *self, avl_node_t *node) 
{
    self->avl_set_node_heap_top++;
    self->avl_set_node_heap[self->avl_set_node_heap_top] = node;
}

static individual_t *
sim_alloc_individual(sim_t *self) 
{
    individual_t *ind;
    if (self->individual_heap_top < 0) {
        fatal_error("Out of individual space");
    }
    ind = self->individual_heap[self->individual_heap_top];
    avl_init_tree(&ind->ancestry, avl_int_map_compare, NULL);
    self->individual_heap_top--;
    return ind;
}

static void
sim_free_individual(sim_t *self, individual_t *ind) 
{
    self->individual_heap_top++;
    self->individual_heap[self->individual_heap_top] = ind;
}

static avl_node_t *
sim_alloc_set_map_node(sim_t *self, unsigned int key) 
{
    avl_node_t *node = NULL;
    set_map_value_t *v;
    if (self->avl_set_map_node_heap_top < 0) {
        fatal_error("Out of set_map space");
    }
    node = self->avl_set_map_node_heap[self->avl_set_map_node_heap_top];
    self->avl_set_map_node_heap_top--;
    v = (set_map_value_t *) node->item;
    v->key = key; 
    avl_init_tree(&v->value, avl_set_compare, NULL);
    return node;
}

static void
sim_free_avl_set_map_node(sim_t *self, avl_node_t *node) 
{
    self->avl_set_map_node_heap_top++;
    self->avl_set_map_node_heap[self->avl_set_map_node_heap_top] = node;
}

static avl_node_t *
sim_alloc_int_map_node(sim_t *self, unsigned int key, unsigned int value) 
{
    avl_node_t *node = NULL;
    int_map_value_t *v;
    if (self->avl_int_map_node_heap_top < 0) {
        fatal_error("Out of int_map space");
    }
    node = self->avl_int_map_node_heap[self->avl_int_map_node_heap_top];
    self->avl_int_map_node_heap_top--;
    v = (int_map_value_t *) node->item;
    v->key = key; 
    v->value = value;
    return node;
}

static void
sim_free_avl_int_map_node(sim_t *self, avl_node_t *node) 
{
    self->avl_int_map_node_heap_top++;
    self->avl_int_map_node_heap[self->avl_int_map_node_heap_top] = node;
}



static void
sim_add_pixel_to_occupancy(sim_t *self, unsigned int h, unsigned int pixel)
{
    avl_node_t *q_node, *h_node, *ret;
    avl_tree_t *tree;
    uint64_t set_value;
    set_map_value_t smv_search;
    set_map_value_t *smv;
    if (h > 0) {
        /* Find the node in Q for occupancy h */
        smv_search.key = h;
        q_node = avl_search(&self->Q, &smv_search);
        if (q_node == NULL) {
            q_node = sim_alloc_set_map_node(self, h);
            avl_insert_node(&self->Q, q_node);
        }
        smv = ((set_map_value_t *) q_node->item);
        tree = &smv->value;
        /* insert pixel into this tree */
        set_value = (uint64_t) pixel;
        h_node = sim_alloc_set_node(self, set_value);
        ret = avl_insert_node(tree, h_node);
        assert(ret != NULL);
    }
}

static void
sim_remove_pixel_from_occupancy(sim_t *self, unsigned int h, unsigned int pixel)
{
    avl_node_t *q_node, *h_node;
    avl_tree_t *tree;
    uint64_t set_value;
    set_map_value_t smv_search;
    set_map_value_t *smv;
    if (h > 0) {
        /* Find the node in Q for occupancy h */
        smv_search.key = h;
        q_node = avl_search(&self->Q, &smv_search);
        assert(q_node != NULL);
        smv = ((set_map_value_t *) q_node->item);
        tree = &smv->value;
        /* remove pixel from this tree */
        set_value = (uint64_t) pixel;
        h_node = avl_search(tree, &set_value);
        assert(h_node != NULL);
        avl_unlink_node(tree, h_node);
        sim_free_avl_set_node(self, h_node); 
        if (avl_count(tree) == 0) {
            avl_unlink_node(&self->Q, q_node);
            sim_free_avl_set_map_node(self, q_node);
        }
    }
}

static void
sim_add_individual_to_pixel(sim_t *self, unsigned int pixel, individual_t *ind)
{
    unsigned int h;
    uintptr_t id = (uintptr_t) ind;
    avl_node_t *node, *ret;
    /* insert the id into this pixel */
    node = sim_alloc_set_node(self, id);
    ret = avl_insert_node(&self->P[pixel], node);
    assert(ret != NULL);
    h = avl_count(&self->P[pixel]);
    sim_add_pixel_to_occupancy(self, h, pixel);
    sim_remove_pixel_from_occupancy(self, h - 1, pixel);
}

static void
sim_add_individual(sim_t *self, individual_t *ind)
{
    unsigned int j, pixel;
    sim_get_disc_pixels(self, ind->location); 
    for (j = 1; j <= self->pixel_buffer[0]; j++) {
        pixel = self->pixel_buffer[j];
        sim_add_individual_to_pixel(self, pixel, ind);
    }
    self->sample_size++;
}


static void 
sim_remove_individual_from_pixel(sim_t *self, unsigned int pixel, 
        individual_t *ind) 
{
    avl_node_t *node;
    uintptr_t id = (uintptr_t) ind;
    unsigned int h;
    /* remove the id from this pixel */
    node = avl_search(&self->P[pixel], &id);
    assert(node != NULL);
    avl_unlink_node(&self->P[pixel], node);
    sim_free_avl_set_node(self, node);
    h = avl_count(&self->P[pixel]);
    sim_add_pixel_to_occupancy(self, h, pixel);
    sim_remove_pixel_from_occupancy(self, h + 1, pixel);
}

static void 
sim_remove_individual(sim_t *self, individual_t *ind)
{
    unsigned int j, pixel;
    sim_get_disc_pixels(self, ind->location); 
    for (j = 1; j <= self->pixel_buffer[0]; j++) {
        pixel = self->pixel_buffer[j];
        sim_remove_individual_from_pixel(self, pixel, ind);
    }
    sim_free_individual(self, ind);
    self->sample_size--;

}
/*
 * Sets up the simulation so that we can enter the main loop.
 */
static void
sim_initialise(sim_t *self)
{
    unsigned int j, k, l, pixel;
    double *x;
    individual_t *ind;
    avl_node_t *node;
    for (j = 0; j < self->n; j++) {
        x = self->X + 2 * j;
        assert(0.0 <= x[0] && x[0] < self->L);
        assert(0.0 <= x[1] && x[1] < self->L);
        ind = sim_alloc_individual(self);
        ind->location[0] = x[0];
        ind->location[1] = x[1];
        for (l = 0; l < self->m; l++) {
            node = sim_alloc_int_map_node(self, l, j + 1);
            avl_insert_node(&ind->ancestry, node);
        }
        sim_get_disc_pixels(self, x); 
        for (k = 1; k <= self->pixel_buffer[0]; k++) {
            pixel = self->pixel_buffer[k];
            sim_add_individual_to_pixel(self, pixel, ind);
        }
    }
    for (l = 0; l < self->m; l++) {
        self->eta[l] = self->n + 1;
        self->coalescence_map[l] = 0;
        for (j = 0; j < 2 * self->n; j++) {
            self->pi[l][j] = 0;
            self->tau[l][j] = 0.0;
        }
    }
    self->sample_size = self->n;
    self->coalesced_loci[0] = 0;
    self->t = 0.0;
    self->ancestral_material = self->n * self->m;
    /* initialise counters */
    self->num_jumps = 0;
    self->num_generated_events = 0;
    for (j = 0; j < gsl_pow_2(self->N); j++) {
        self->total_E_size[j] = 0;
        self->total_H_size[j] = 0;
        self->total_S_size[j] = 0;
        self->total_C_size[j] = 0;
        self->total_G_size[j] = 0;
        
    }
}
static double 
sim_ubar(sim_t *self, unsigned int n) 
{
    assert(n <= self->max_occupancy);
    return self->ubar[n];
}

/* 
 * Selects a parent for the next locus which is gap loci distant.
 */
static unsigned int 
sim_select_parent(sim_t *self, unsigned int current_parent, unsigned int gap)
{
    unsigned int ret = current_parent;
    double nu = self->nu;
    unsigned int k = gap;
    double rho = self->rho;
    double psi = ((nu - 1) / nu) 
            * (1 - gsl_pow_int(1 - nu * rho / (nu - 1), k));
    if (gsl_rng_uniform(self->rng) < psi) {
        if (nu == 2) {
            ret = (ret + 1) % 2;
        } else {
            ret = current_parent;
            while (ret == current_parent) {
                ret = gsl_rng_uniform_int(self->rng, nu);
            }
        }

    }
    return ret;
}

static void 
sim_update_parental_ancestry(sim_t *self, unsigned int parent, 
        unsigned int locus, unsigned int ancestry)
{
    avl_node_t *node, *ret;
    individual_t *ind;
    int_map_value_t *imv, imv_search;
    ind = self->parent_buffer[parent];
    node = sim_alloc_int_map_node(self, locus, ancestry);
    ret = avl_insert_node(&ind->ancestry, node);
    if (ret == NULL) {
        /* Coalescence */
        sim_free_avl_int_map_node(self, node);
        /* Get the existing node */
        imv_search.key = locus;
        node = avl_search(&ind->ancestry, &imv_search);
        assert(node != NULL);
        imv = (int_map_value_t *) node->item;
        if (self->coalescence_map[locus] == 0) {
            self->coalescence_map[locus] = 1;
            self->pi[locus][imv->value] = self->eta[locus];
            self->tau[locus][self->eta[locus]] = self->t;
            imv->value = self->eta[locus];
            self->eta[locus]++;
            self->coalesced_loci[0]++;
            self->coalesced_loci[self->coalesced_loci[0]] = locus;
        }
        self->pi[locus][ancestry] = self->eta[locus] - 1;
        self->ancestral_material--;
    }
}

/* 
 * Generates parental individuals and inserts them into the parent
 * buffer based on the individuals in the child_buffer. 
 */
static void
sim_generate_parents(sim_t *self, unsigned int C_size)
{
    unsigned int j, parent, previous_locus;
    individual_t **C = self->child_buffer;
    avl_node_t *node;
    int_map_value_t *locus_mapping;
    for (j = 0; j < self->nu; j++) {
        self->parent_buffer[j] = sim_alloc_individual(self);
    }
    for (j = 0; j < C_size; j++) {
        /* Assign the parent for the first locus */
        parent = gsl_rng_uniform_int(self->rng, self->nu);
        node = C[j]->ancestry.head;
        assert(node != NULL);
        locus_mapping = (int_map_value_t *) node->item;
        sim_update_parental_ancestry(self, parent, locus_mapping->key,
                locus_mapping->value);
        previous_locus = locus_mapping->key;
        while (node->next != NULL) {
            node = node->next; 
            locus_mapping = (int_map_value_t *) node->item;
            parent = sim_select_parent(self, parent, locus_mapping->key
                    - previous_locus);
            sim_update_parental_ancestry(self, parent, locus_mapping->key,
                    locus_mapping->value);
            previous_locus = locus_mapping->key;
        }
        /* Now free the ancestry nodes */
        for (node = C[j]->ancestry.head; node != NULL; node = node->next) {
            sim_free_avl_int_map_node(self, node);
        }
    }
        
}

static void 
sim_run(sim_t *self)
{
    unsigned int j, k, pixel, max_occupancy, occupancy, pixel_count, v[2];
    double Lambda_const = self->lambda * gsl_pow_2(self->s / self->L);
    double z[] = {0.0, 0.0};
    double *p = self->probability_buffer;
    individual_t **S = self->intersected_buffer;
    individual_t **C = self->child_buffer;
    unsigned int S_size, C_size;
    uint64_t set_item;
    double Lambda, jump_proba, *x;
    avl_node_t *node;
    set_map_value_t *smv, smv_search;
    individual_t *ind;
    while (self->ancestral_material > self->m && self->t < self->max_time
            && self->num_jumps < self->max_jumps) {
        //sim_print_state(self, 0);
        Lambda = 0.0;
        max_occupancy = 0;
        memset(p, 0, self->max_occupancy * sizeof(double)); 
        for (node = self->Q.head; node != NULL; node = node->next) {
            smv = (set_map_value_t *) node->item;
            occupancy = smv->key;
            pixel_count = avl_count(&smv->value);
            assert(occupancy <= self->max_occupancy);
            p[occupancy - 1] = pixel_count * sim_ubar(self, occupancy);
            Lambda += p[occupancy - 1]; 
            if (occupancy > max_occupancy) {
                max_occupancy = occupancy;
            }
        }
        for (j = 1; j <= max_occupancy; j++) {
            p[j - 1] /= Lambda;
        }
        do {
            self->t += gsl_ran_exponential(self->rng, 
                    1.0 / (Lambda * Lambda_const));
            occupancy = 1 + probability_list_select(p, max_occupancy, 
                    gsl_rng_uniform(self->rng));
            /* choose a pixel with this occupancy uniformly */
            smv_search.key = occupancy;
            node = avl_search(&self->Q, &smv_search);
            assert(node != NULL);
            smv = (set_map_value_t *) node->item;
            k = gsl_rng_uniform_int(self->rng, avl_count(&smv->value));
            node = avl_at(&smv->value, k);
            assert(node != NULL);
            set_item =  *((uint64_t *) node->item);
            pixel = (unsigned int) set_item; 
            index_to_pixel_coord(pixel, self->N, v);
            z[0] = self->s * (v[0] + gsl_rng_uniform(self->rng));
            z[1] = self->s * (v[1] + gsl_rng_uniform(self->rng));
            assert(0.0 <= z[0] && self->L > z[0]);
            assert(0.0 <= z[1] && self->L > z[1]);
            S_size = 0;
            for (node = self->P[pixel].head; node != NULL; node = node->next) {
                set_item =  *((uint64_t *) node->item);
                ind = (individual_t *) set_item; 
                x = ind->location; 
                //printf("\tchecking %p (%f, %f): %f\n", ind, x[0], x[1], 
                //        sqrt(torus_squared_distance(z, x, self->L)));
                if (torus_squared_distance(z, x, self->L) < gsl_pow_2(self->r)) {
                    S[S_size] = ind;
                    S_size++;
                }
            }
            jump_proba = sim_ubar(self, S_size) / sim_ubar(self, occupancy); 
            //printf("|S| = %d: o = %d; jump proba = %f\n", S_size, occupancy, jump_proba);
            assert(jump_proba <= 1.0);
        } while (gsl_rng_uniform(self->rng) < 1.0 - jump_proba);
        C_size = gsl_ran_discrete(self->rng, self->beta_distributions[S_size]);
        gsl_ran_choose(self->rng, C, C_size, S, S_size, 
                sizeof(individual_t *));
        /*
        printf("Children: %d\n", C[0]);
        for (j = 1; j <= C[0]; j++) {
            sim_print_individual(self, C[j]);
        }
        */
        sim_generate_parents(self, C_size);
        for (j = 0; j < C_size; j++) {
            sim_remove_individual(self, C[j]);
        }
        /* add in the parents, if they have any ancestral material */
        //printf("Parents:\n");
        for (j = 0; j < self->nu; j++) {
            ind = self->parent_buffer[j];
            if (avl_count(&ind->ancestry) == 0) {
                sim_free_individual(self, ind);
            } else {
                self->total_G_size[pixel]++;
                random_point_torus_disc(ind->location, z, self->r, self->L, 
                        self->rng);
                sim_add_individual(self, ind);
                //sim_print_individual(self, self->next_id - 1); 
            }
        }

        /* Clear the coalescence markers */
        /*
        if (self->coalesced_loci[0] > 0) {
            sim_print_state(self);
             
        }
        */
        for (j = 1; j <= self->coalesced_loci[0]; j++) {
            k = self->coalesced_loci[j];
            self->coalescence_map[k] = 0;
        }
        self->coalesced_loci[0] = 0;

        //printf("%f\n", self->t);
    } 
}


