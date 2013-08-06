/*
** Copyright (C) 2013 Jerome Kelleher <jerome.kelleher@ed.ac.uk>
**  
** This file is part of discsim.
** 
** discsim is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
** 
** discsim is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with discsim.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "sim.h"
#include "util.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>


const char * 
sim_error_message(int err)
{
    const char *ret = "Unknown error";
    switch (err) {
        case ERR_ALLOC_FAILED:
            ret = "Memory allocation failure";
            break;
        case ERR_BAD_PIXEL_SIZE:
            ret = "Bad pixel size";
            break;
        case ERR_OUT_OF_INT_MAP_NODES:
            ret = "Out of int map nodes";
            break;
        case ERR_OUT_OF_SET_MAP_NODES:
            ret = "Out of set map nodes";
            break;
        case ERR_OUT_OF_AVL_SET_NODES:
            ret = "Out of avl set nodes";
        case ERR_OUT_OF_INDIVIDUALS:
            ret = "Out of individuals";
            break;
        case ERR_AVL_OP_FAILED:
            ret = "AVL tree operation failed";
            break;
    }
    return ret;
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
 * with num_parents parents and impact u.
 */
double 
equilibrium_density(double num_parents, double u)
{
    double N = num_parents / u;
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
torus_squared_distance(const double *p1, const double *p2, double torus_diameter) 
{
    double xabs = fabs(p2[0] - p1[0]);
    double yabs = fabs(p2[1] - p1[1]);
    double xd = GSL_MIN(xabs, torus_diameter - xabs);
    double yd = GSL_MIN(yabs, torus_diameter - yabs);
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
 * side torus_diameter.
 *
 * This uses the canonical method for generating a point uniformly distributed within
 * the unit circle (from TAOCP), then scales and takes the modulus appropriately.
 */
static inline void 
random_point_torus_disc(double *p, const double *centre, const double r, 
        const double torus_diameter, gsl_rng *generator)
{
    register double pixel_size = 1.1;
    register double x;
    register double y;
    while (pixel_size >= 1) {   
        x = 2 * gsl_rng_uniform(generator) - 1;
        y = 2 * gsl_rng_uniform(generator) - 1;
        pixel_size = (x * x) + (y * y);
    }
    x = centre[0] + (x * r);
    y = centre[1] + (y * r);
    p[0] = fmod(x + torus_diameter, torus_diameter);
    p[1] = fmod(y + torus_diameter, torus_diameter); 
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


/* Returns 1 if a disc of radius r intersects with a pixel of side pixel_size
 * with bottom-left corner x*pixel_size, y*pixel_size.
 */
static int 
disc_intersects_pixel(double *z, int x, int y, double r, double pixel_size)
{
    // Get the closest point in the pixel to the centre of the disc
    double closest_x = clamp(z[0], x * pixel_size, (x + 1) * pixel_size);
    double closest_y = clamp(z[1], y * pixel_size, (y + 1) * pixel_size);
    // Calculate the distance between the disc'pixel_size center and this closest point
    double dx2 = gsl_pow_2(z[0] - closest_x);
    double dy2 = gsl_pow_2(z[1] - closest_y);
    return  dx2 + dy2 < gsl_pow_2(r);
}


/* Gets the pixels for the disc at the specified location and updates the pixel 
 * buffer.
 */
static void 
get_pixels_general(double r, double pixel_size, int N, double *z, unsigned int *p)
{
    int x_min = (int) floor((z[0] - r) / pixel_size);
    int x_max = (int) ceil((z[0] + r) / pixel_size);
    int y_min = (int) floor((z[1] - r) / pixel_size);
    int y_max = (int) ceil((z[1] + r) / pixel_size);
    int x, y;
    int coord[2];
    p[0] = 0;
    for (y = y_min; y <= y_max; y++) {
        for (x = x_min; x <= x_max; x++) {
            if (disc_intersects_pixel(z, x, y, r, pixel_size)) {
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
get_pixels(double r, double pixel_size, int N, double *z, unsigned int *p)
{
    double r2 = r * r;
    double a[2], b[2], c[2], d[2], x[2];
    double torus_diameter = pixel_size * N;
    int v[2];
    assert(pixel_size >= 1.0);
    v[0] = (unsigned int) floor(z[0] / pixel_size); 
    v[1] = (unsigned int) floor(z[1] / pixel_size); 
    p[0] = 1;
    p[p[0]] = pixel_coord_to_index(v[0], v[1], N);
    a[0] = pixel_size * v[0];
    a[1] = pixel_size * v[1];
    b[0] = a[0] + r;
    b[1] = a[1] + r;
    d[0] = a[0] + pixel_size;
    d[1] = a[1] + pixel_size;
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
        if (torus_squared_distance(a, z, torus_diameter) < r2) {
            p[0]++;
            p[p[0]] = pixel_coord_to_index((v[0] - 1 + N ) % N, 
                    (v[1] - 1 + N) % N, N);
        } 
    } 
    if (z[0] < b[0] && z[1] > c[1]) {
        x[0] = a[0];
        x[1] = d[1];
        if (torus_squared_distance(x, z, torus_diameter) < r2) {
            p[0]++;
            p[p[0]] = pixel_coord_to_index((v[0] - 1 + N) % N, 
                    (v[1] + 1) % N, N);
        }
    } 
    if (z[0] > c[0] && z[1] > c[1]) {
        if (torus_squared_distance(d, z, torus_diameter) < r2) {
            p[0]++;
            p[p[0]] = pixel_coord_to_index((v[0] + 1) % N, 
                    (v[1] + 1) % N, N);
        }
    } 
    if (z[0] > c[0] && z[1] < b[1]) {
        x[0] = d[0];
        x[1] = a[1];
        if (torus_squared_distance(x, z, torus_diameter) < r2) {
            p[0]++;
            p[p[0]] = pixel_coord_to_index((v[0] + 1) % N, 
                    (v[1] - 1 + N) % N, N);
        }
    }
#ifdef CHECK_PIXELS 
    unsigned int j;
    unsigned int p2[10];
    get_pixels_general(r, pixel_size, N, z, p2);
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
    double r = self->event_classes[0].r;
    if (self->pixel_size < 1) {
        get_pixels_general(r, self->pixel_size, self->N, z, self->pixel_buffer);
    } else {
        get_pixels(r, self->pixel_size, self->N, z, self->pixel_buffer);
    }

}


/*
 * Sets the max occupancy to reflect the parameters of the simulation 
 * using the equilibrium prediction, using the specified fraction of headroom
 * over this estimate.
 */
void
sim_set_max_occupancy(sim_t *self, double headroom)
{
    double u = self->event_classes[0].u;
    double delta = equilibrium_density(self->num_parents, u);
    double o = delta * (gsl_pow_2(self->pixel_size) + 4 * self->pixel_size + M_PI);
    o += o * headroom;
    self->max_occupancy = (unsigned int) o; 

}

int
sim_alloc(sim_t *self)
{
    int ret = 0;
    unsigned int j, k, max;
    double *b = NULL;
    double r = self->event_classes[0].r;
    double u = self->event_classes[0].u;
    self->N = self->torus_diameter / self->pixel_size;
    if (fmod(self->torus_diameter, self->pixel_size) != 0.0) {
        ret = ERR_BAD_PIXEL_SIZE;
        goto out;
    }
    self->rng = gsl_rng_alloc(gsl_rng_default);
    if (self->rng == NULL) {
        ret = ERR_ALLOC_FAILED;
        goto out;
    }
    gsl_rng_set(self->rng, self->random_seed);
    self->P = xmalloc(gsl_pow_2(self->N) * sizeof(avl_tree_t));
    /* set up buffers */
    self->max_disc_pixels = gsl_pow_2(4 * ((int)(r / self->pixel_size) + 1));
    self->pixel_buffer = xmalloc(self->max_disc_pixels * sizeof(unsigned int));
    self->probability_buffer = xmalloc(self->max_occupancy * sizeof(double));
    self->intersected_buffer = xmalloc(
            (self->max_occupancy + 1) * sizeof(unsigned int));
    self->child_buffer = xmalloc(
            (self->max_occupancy + 1) * sizeof(unsigned int));
    self->parent_buffer = xmalloc((self->num_parents + 1) * sizeof(avl_node_t *));
    /* Set up the AVL trees */
    avl_init_tree(&self->Q, avl_set_map_compare, NULL);
    for (j = 0; j < gsl_pow_2(self->N); j++) {
        avl_init_tree(&self->P[j], avl_set_compare, NULL);
        
    }   
    /* Set up the set memory */
    max = self->max_disc_pixels * self->max_population_size;
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
    max = self->max_population_size;
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
    max = (self->sample_size * gsl_max(self->num_parents, 2)) * self->num_loci;
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
    self->eta = xmalloc(self->num_loci * sizeof(unsigned int));
    self->pi = xmalloc(self->num_loci * sizeof(int *));
    self->tau = xmalloc(self->num_loci * sizeof(double *));
    self->pi_mem = xmalloc(2 * self->sample_size * self->num_loci * sizeof(int));
    self->tau_mem = xmalloc(2 * self->sample_size * self->num_loci * sizeof(double));
    for (j = 0; j < self->num_loci; j++) {
        self->pi[j] = self->pi_mem + j * 2 * self->sample_size; 
        self->tau[j] = self->tau_mem + j * 2 * self->sample_size; 
    }
    self->coalescence_map = xmalloc(self->num_loci * sizeof(unsigned int));
    self->coalesced_loci = xmalloc((self->num_loci + 1) * sizeof(unsigned int));
    /* Precalculate ubar */
    self->ubar = xmalloc((self->max_occupancy + 1) * sizeof(double));
    for (j = 0; j <= self->max_occupancy; j++) {
        self->ubar[j] = 1.0 - gsl_pow_int(1.0 - u, j);
    }
    /* Precalculate beta */
    b = xmalloc((self->max_occupancy + 1) * sizeof(double));
    b[0] = 0.0;
    self->beta_distributions = xmalloc((self->max_occupancy + 1) 
            * sizeof(gsl_ran_discrete_t *));
    for (k = 1; k <= self->max_occupancy; k++) {
        /* first figure out how many places we need to get below threshold */
        j = 1;
        while (beta(k, j, u) > self->beta_threshold && j < k) {
            j++;
        }
        max = j;
        for (j = 1; j <= max; j++) {
            b[j] = beta(k, j, u);
        }
        self->beta_distributions[k] = gsl_ran_discrete_preproc(max + 1, b);
        if (self->beta_distributions[k] == NULL) {
            ret = ERR_ALLOC_FAILED;
            goto out;
        }
    }
out:
    if (b != NULL) {
        free(b);
    }
    return ret;
}

void
sim_free(sim_t *self)
{
    unsigned int j;
    if (self->rng != NULL) {
        gsl_rng_free(self->rng);
    }
    xfree(self->P);
    xfree(self->pixel_buffer);
    xfree(self->probability_buffer);
    xfree(self->intersected_buffer);
    xfree(self->child_buffer);
    xfree(self->parent_buffer);
    xfree(self->avl_set_value_mem);
    xfree(self->avl_set_node_mem);
    xfree(self->avl_set_node_heap);
    xfree(self->individual_mem);
    xfree(self->individual_heap);
    xfree(self->avl_set_map_value_mem);
    xfree(self->avl_set_map_node_mem);
    xfree(self->avl_set_map_node_heap);
    xfree(self->avl_int_map_value_mem);
    xfree(self->avl_int_map_node_mem);
    xfree(self->avl_int_map_node_heap);
    xfree(self->pi);
    xfree(self->tau);
    xfree(self->eta);
    xfree(self->pi_mem);
    xfree(self->tau_mem);
    xfree(self->coalescence_map);
    xfree(self->coalesced_loci);
    xfree(self->ubar);
    if (self->beta_distributions != NULL) {
        for (j = 1; j < self->max_occupancy + 1; j++) {
            if (self->beta_distributions[j] != NULL) {
                gsl_ran_discrete_free(self->beta_distributions[j]);
            }
        }
        free(self->beta_distributions);
    }
}

static avl_node_t *
sim_alloc_avl_set_node(sim_t *self, uint64_t value) 
{
    avl_node_t *node = NULL;
    if (self->avl_set_node_heap_top < 0) {
        goto out; 
    }
    node = self->avl_set_node_heap[self->avl_set_node_heap_top];
    self->avl_set_node_heap_top--;
    *((uint64_t *) node->item) = value; 
out:
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
    individual_t *ind = NULL;
    if (self->individual_heap_top < 0) {
        goto out; 
    }
    ind = self->individual_heap[self->individual_heap_top];
    avl_init_tree(&ind->ancestry, avl_int_map_compare, NULL);
    self->individual_heap_top--;
out:
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
         goto out;
    }
    node = self->avl_set_map_node_heap[self->avl_set_map_node_heap_top];
    self->avl_set_map_node_heap_top--;
    v = (set_map_value_t *) node->item;
    v->key = key; 
    avl_init_tree(&v->value, avl_set_compare, NULL);
out:
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
        goto out; 
    }
    node = self->avl_int_map_node_heap[self->avl_int_map_node_heap_top];
    self->avl_int_map_node_heap_top--;
    v = (int_map_value_t *) node->item;
    v->key = key; 
    v->value = value;
out:
    return node;
}

static void
sim_free_avl_int_map_node(sim_t *self, avl_node_t *node) 
{
    self->avl_int_map_node_heap_top++;
    self->avl_int_map_node_heap[self->avl_int_map_node_heap_top] = node;
}



static int 
sim_add_pixel_to_occupancy(sim_t *self, unsigned int h, unsigned int pixel)
{
    int ret = 0;
    avl_node_t *q_node, *h_node;
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
            if (q_node == NULL) {
                ret = ERR_OUT_OF_SET_MAP_NODES;
                goto out;
            }
            avl_insert_node(&self->Q, q_node);
        }
        smv = ((set_map_value_t *) q_node->item);
        tree = &smv->value;
        /* insert pixel into this tree */
        set_value = (uint64_t) pixel;
        h_node = sim_alloc_avl_set_node(self, set_value);
        if (h_node == NULL) {
            ret = ERR_OUT_OF_AVL_SET_NODES;
            goto out;
        }
        if (avl_insert_node(tree, h_node) == NULL) {
            ret = ERR_AVL_OP_FAILED;
            goto out;
        }
    }
out:
    return ret;
}

static int 
sim_remove_pixel_from_occupancy(sim_t *self, unsigned int h, unsigned int pixel)
{
    int ret = 0;
    avl_node_t *q_node, *h_node;
    avl_tree_t *tree;
    uint64_t set_value;
    set_map_value_t smv_search;
    set_map_value_t *smv;
    if (h > 0) {
        /* Find the node in Q for occupancy h */
        smv_search.key = h;
        q_node = avl_search(&self->Q, &smv_search);
        if (q_node == NULL) {
            ret = ERR_AVL_OP_FAILED;
            goto out;
        }
        smv = ((set_map_value_t *) q_node->item);
        tree = &smv->value;
        /* remove pixel from this tree */
        set_value = (uint64_t) pixel;
        h_node = avl_search(tree, &set_value);
        if (h_node == NULL) {
            ret = ERR_AVL_OP_FAILED;
            goto out;
        }
        avl_unlink_node(tree, h_node);
        sim_free_avl_set_node(self, h_node); 
        if (avl_count(tree) == 0) {
            avl_unlink_node(&self->Q, q_node);
            sim_free_avl_set_map_node(self, q_node);
        }
    }
out:
    return ret;
}

static int 
sim_add_individual_to_pixel(sim_t *self, unsigned int pixel, individual_t *ind)
{
    int ret = 0;
    unsigned int h;
    uintptr_t id = (uintptr_t) ind;
    avl_node_t *node;
    /* insert the id into this pixel */
    node = sim_alloc_avl_set_node(self, id);
    if (avl_insert_node(&self->P[pixel], node) == NULL) {
        ret = ERR_AVL_OP_FAILED;
        goto out;
    }
    h = avl_count(&self->P[pixel]);
    ret = sim_add_pixel_to_occupancy(self, h, pixel);
    if (ret != 0) {
        goto out;
    }
    ret = sim_remove_pixel_from_occupancy(self, h - 1, pixel);
out:
    return ret;
}

static int 
sim_add_individual(sim_t *self, individual_t *ind)
{
    int ret = 0;
    unsigned int j, pixel;
    sim_get_disc_pixels(self, ind->location); 
    for (j = 1; j <= self->pixel_buffer[0]; j++) {
        pixel = self->pixel_buffer[j];
        ret = sim_add_individual_to_pixel(self, pixel, ind);
        if (ret != 0) {
            goto out;
        }
    }
    self->sample_size++;
out:
    return ret;
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
int
sim_initialise(sim_t *self)
{
    int ret = 0;
    unsigned int j, k, l, pixel;
    double *x;
    individual_t *ind;
    avl_node_t *node;
    for (j = 0; j < self->sample_size; j++) {
        x = self->sample + 2 * j;
        assert(0.0 <= x[0] && x[0] < self->torus_diameter);
        assert(0.0 <= x[1] && x[1] < self->torus_diameter);
        ind = sim_alloc_individual(self);
        if (ind == NULL) {
            ret = ERR_OUT_OF_INDIVIDUALS;
            goto out;
        }
        ind->location[0] = x[0];
        ind->location[1] = x[1];
        for (l = 0; l < self->num_loci; l++) {
            node = sim_alloc_int_map_node(self, l, j + 1);
            if (node == NULL) {
                ret = ERR_OUT_OF_INT_MAP_NODES;
                goto out;
            }
            if (avl_insert_node(&ind->ancestry, node) == NULL) {
                ret = ERR_AVL_OP_FAILED;
                goto out;
            }
        }
        sim_get_disc_pixels(self, x); 
        for (k = 1; k <= self->pixel_buffer[0]; k++) {
            pixel = self->pixel_buffer[k];
            ret = sim_add_individual_to_pixel(self, pixel, ind);
            if (ret != 0) {
                goto out;
            }
        }
    }
    for (l = 0; l < self->num_loci; l++) {
        self->eta[l] = self->sample_size + 1;
        self->coalescence_map[l] = 0;
        for (j = 0; j < 2 * self->sample_size; j++) {
            self->pi[l][j] = 0;
            self->tau[l][j] = 0.0;
        }
    }
    self->sample_size = self->sample_size;
    self->coalesced_loci[0] = 0;
    self->time = 0.0;
    self->successful_events = 0;
    self->ancestral_material = self->sample_size * self->num_loci;
out:
    return ret;
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
    double num_parents = self->num_parents;
    unsigned int k = gap;
    double rho = self->recombination_probability;
    double psi = ((num_parents - 1) / num_parents) 
            * (1 - gsl_pow_int(1 - num_parents * rho / (num_parents - 1), k));
    if (gsl_rng_uniform(self->rng) < psi) {
        if (num_parents == 2) {
            ret = (ret + 1) % 2;
        } else {
            ret = current_parent;
            while (ret == current_parent) {
                ret = gsl_rng_uniform_int(self->rng, num_parents);
            }
        }

    }
    return ret;
}

static int  
sim_update_parental_ancestry(sim_t *self, unsigned int parent, 
        unsigned int locus, unsigned int ancestry)
{
    int ret = 0;
    avl_node_t *node, *avl_ret;
    individual_t *ind;
    int_map_value_t *imv, imv_search;
    ind = self->parent_buffer[parent];
    node = sim_alloc_int_map_node(self, locus, ancestry);
    if (node == NULL) {
        ret = ERR_OUT_OF_INT_MAP_NODES;
        goto out;
    }
    avl_ret = avl_insert_node(&ind->ancestry, node);
    if (avl_ret == NULL) {
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
            self->tau[locus][self->eta[locus]] = self->time;
            imv->value = self->eta[locus];
            self->eta[locus]++;
            self->coalesced_loci[0]++;
            self->coalesced_loci[self->coalesced_loci[0]] = locus;
        }
        self->pi[locus][ancestry] = self->eta[locus] - 1;
        self->ancestral_material--;
    }
out:
    return ret;
}

/* 
 * Generates parental individuals and inserts them into the parent
 * buffer based on the individuals in the child_buffer. 
 */
static int 
sim_generate_parents(sim_t *self, unsigned int C_size)
{
    int ret = 0;
    unsigned int j, parent, previous_locus;
    individual_t **C = self->child_buffer;
    avl_node_t *node;
    int_map_value_t *locus_mapping;
    for (j = 0; j < self->num_parents; j++) {
        self->parent_buffer[j] = sim_alloc_individual(self);
        if (self->parent_buffer[j] == NULL) {
            ret = ERR_OUT_OF_INDIVIDUALS;
            goto out;
        }
    }
    for (j = 0; j < C_size; j++) {
        /* Assign the parent for the first locus */
        parent = gsl_rng_uniform_int(self->rng, self->num_parents);
        node = C[j]->ancestry.head;
        assert(node != NULL);
        locus_mapping = (int_map_value_t *) node->item;
        ret = sim_update_parental_ancestry(self, parent, locus_mapping->key,
                locus_mapping->value);
        if (ret != 0) {
            goto out;
        }
        previous_locus = locus_mapping->key;
        while (node->next != NULL) {
            node = node->next; 
            locus_mapping = (int_map_value_t *) node->item;
            parent = sim_select_parent(self, parent, locus_mapping->key
                    - previous_locus);
            ret = sim_update_parental_ancestry(self, parent, locus_mapping->key,
                    locus_mapping->value);
            if (ret != 0) {
                goto out;
            }
            previous_locus = locus_mapping->key;
        }
        /* Now free the ancestry nodes */
        for (node = C[j]->ancestry.head; node != NULL; node = node->next) {
            sim_free_avl_int_map_node(self, node);
        }
    }
out:
    return ret;
}

int 
sim_simulate(sim_t *self, unsigned int max_jumps, double max_time)
{
    int ret = 0;
    unsigned int j, k, pixel, max_occupancy, occupancy, pixel_count, v[2];
    double r = self->event_classes[0].rate;
    double lambda = self->event_classes[0].rate;
    double Lambda_const = lambda * gsl_pow_2(self->pixel_size / self->torus_diameter);
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
    while (self->ancestral_material > self->num_loci && self->time < max_time
            && self->successful_events < max_jumps) {
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
            self->time += gsl_ran_exponential(self->rng, 
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
            z[0] = self->pixel_size * (v[0] + gsl_rng_uniform(self->rng));
            z[1] = self->pixel_size * (v[1] + gsl_rng_uniform(self->rng));
            assert(0.0 <= z[0] && self->torus_diameter > z[0]);
            assert(0.0 <= z[1] && self->torus_diameter > z[1]);
            S_size = 0;
            for (node = self->P[pixel].head; node != NULL; node = node->next) {
                set_item =  *((uint64_t *) node->item);
                ind = (individual_t *) set_item; 
                x = ind->location; 
                if (torus_squared_distance(z, x, self->torus_diameter) < gsl_pow_2(r)) {
                    S[S_size] = ind;
                    S_size++;
                }
            }
            jump_proba = sim_ubar(self, S_size) / sim_ubar(self, occupancy); 
            assert(jump_proba <= 1.0);
        } while (gsl_rng_uniform(self->rng) < 1.0 - jump_proba);
        self->successful_events++;
        C_size = gsl_ran_discrete(self->rng, self->beta_distributions[S_size]);
        gsl_ran_choose(self->rng, C, C_size, S, S_size, sizeof(individual_t *));
        ret = sim_generate_parents(self, C_size);
        if (ret != 0) {
            goto out;
        }
        for (j = 0; j < C_size; j++) {
            sim_remove_individual(self, C[j]);
        }
        /* add in the parents, if they have any ancestral material */
        for (j = 0; j < self->num_parents; j++) {
            ind = self->parent_buffer[j];
            if (avl_count(&ind->ancestry) == 0) {
                sim_free_individual(self, ind);
            } else {
                random_point_torus_disc(ind->location, z, r, self->torus_diameter, 
                        self->rng);
                ret = sim_add_individual(self, ind);
                if (ret != 0) {
                    goto out;
                }
            }
        }
        /* Clear the coalescence markers */
        for (j = 1; j <= self->coalesced_loci[0]; j++) {
            k = self->coalesced_loci[j];
            self->coalescence_map[k] = 0;
        }
        self->coalesced_loci[0] = 0;
    }
out:
    return ret;
}


void
sim_print_parameters(sim_t *self)
{
    unsigned int j;
    double *x;
    printf("# torus_diameter = %f\n", self->torus_diameter);
    printf("# num_loci = %u\n", self->num_loci);
    printf("# pixel_size = %f\n", self->pixel_size);
    printf("# num_parents = %u\n", self->num_parents);
    printf("# recombination_probability = %f\n", self->recombination_probability);
    printf("# events: \n");
    for (j = 0; j < self->num_event_classes; j++) {
        printf("# event class %d:\n", j);
        printf("# \tu = %f\n", self->event_classes[j].u);
        printf("# \tr = %f\n", self->event_classes[j].r);
        printf("# \trate = %f\n", self->event_classes[j].rate);
    }
    printf("# sample = ");
    for (j = 0; j < self->sample_size; j++) {
        x = self->sample + 2 * j;
        printf("(%f, %f) ", x[0], x[1]); 
    }
    printf("\n# max_population_size = %u\n", self->max_population_size);
    printf("# max_occupancy = %u\n", self->max_occupancy);
    printf("# random_seed = %lu\n", self->random_seed);

}

void
sim_print_state(sim_t *self, int detail)
{
    unsigned int j, k, l, h, v[2];
    avl_tree_t chi;
    double *x;
    uint64_t id, pixel;
    set_map_value_t *set_map_value, search;
    avl_node_t *node, *chi_node;
    avl_init_tree(&chi, avl_set_compare, NULL); 
    individual_t *ind;
    node = NULL; 
    /* fill chi first */
    for (pixel = 0; pixel < gsl_pow_2(self->N); pixel++) {
        for (node = self->P[pixel].head; node != NULL; node = node->next) {
            id = *((uint64_t *) node->item);
            chi_node = sim_alloc_avl_set_node(self, id); 
            if (avl_insert_node(&chi, chi_node) == NULL) {
                sim_free_avl_set_node(self, chi_node);
            }
        }
    }
    assert(avl_count(&chi) == self->sample_size);
    printf("chi = \n");
    for (chi_node = chi.head; chi_node != NULL; chi_node = chi_node->next) {
        id = *((uint64_t *) chi_node->item);
        ind = (individual_t *) id;
        x = ind->location;
        printf("\t %p -> (%f, %f): %d\n", ind, x[0], x[1], 
                avl_count(&ind->ancestry));                
    }
    if (detail > 0) {
        printf("    ");
        for (k = 0; k < self->N; k++) {
            printf("%02d ", k);
        }
        printf("\n    ");
        for (k = 0; k < self->N; k++) {
            printf("___");
        }
        printf("\n");
        for (j = 0; j < self->N; j++) {
            printf(" %02d|", j);
            for (k = 0; k < self->N; k++) {
                v[0] = j;
                v[1] = k;
                pixel = pixel_coord_to_index(v[0], v[1], self->N); 
                h = avl_count(&self->P[pixel]);
                if (h == 0) {
                    printf("   ");
                } else {
                    printf("%02d ", h);
                    /*  Verify the pixel is in Q[h] */
                    search.key = h;
                    node = avl_search(&self->Q, &search);
                    assert(node != NULL);
                    set_map_value = (set_map_value_t *) node->item;
                    node = avl_search(&set_map_value->value, &pixel);
                    assert(node != NULL);
                }
            }
            printf("|\n");
        }
        printf("    ");
        for (k = 0; k < self->N; k++) {
            printf("___");
        }
        printf("\n");
        printf("Q = \n");
        for (j = 0; j < avl_count(&self->Q); j++) {
            node = avl_at(&self->Q, j);
            set_map_value = (set_map_value_t *) node->item;
            h = set_map_value->key;
            printf("\t %d -> ", h);
            for (k = 0; k < avl_count(&set_map_value->value); k++) {
                node = avl_at(&set_map_value->value, k);
                pixel = *((uint64_t *) node->item);
                printf(" %ld, ", pixel);
                /* verify that the occupancy of this pixel is h */
                assert(avl_count(&self->P[pixel]) == h);
            }
            printf("\n");
        }
        printf("Memory:\n");
        printf("\tavl_set_node_heap_top =            %d\n", 
                self->avl_set_node_heap_top);
        printf("\tindividual_heap_top = %d\n", 
                self->individual_heap_top);
        printf("\tavl_set_map_node_heap_top =        %d\n", 
                self->avl_set_map_node_heap_top);
        printf("\tavl_int_map_node_heap_top =        %d\n", 
                self->avl_int_map_node_heap_top);
    }
    if (detail > 1) {
        printf("Ancestry:\n");
        for (l = 0; l < self->num_loci; l++) {
            printf("\t%d (%d)\tpi = ", l, self->eta[l]);
            for (j = 0; j < 2 * self->sample_size; j++) {
                printf("%d ", self->pi[l][j]);
            }
            printf("\ttau = ");
            for (j = 0; j < 2 * self->sample_size; j++) {
                printf("%f ", self->tau[l][j]);
            }
            printf("\n");
        }
    }
    for (chi_node = chi.head; chi_node != NULL; chi_node = chi_node->next) {
        sim_free_avl_set_node(self, chi_node);
    }
}

