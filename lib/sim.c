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

#include <gsl/gsl_sort.h>

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

/* Gets the pixels for the disc at the specified location and updates the pixel 
 * buffer.
 */
static void 
sim_get_disc_pixels(sim_t *self, double *z, double r)
{
    get_pixels_general(r, self->pixel_size, self->N, z, self->pixel_buffer);
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
    double u = self->event_classes[0].u;
    double r;
    if (self->dimension < 1 || self->dimension > 2) {
        ret = ERR_BAD_DIMENSION;
        goto out;
    }
    if (self->dimension == 1 && self->pixel_size != 2.0) {
        ret = ERR_BAD_PIXEL_SIZE;
        goto out;
    }
    self->N = (unsigned int) self->torus_diameter / self->pixel_size;
    /* we assume that uint64_t >= size of a pointer */
    assert(sizeof(uint64_t) >= sizeof(uintptr_t));
    if (fmod(self->torus_diameter, self->pixel_size) != 0.0) {
        ret = ERR_BAD_PIXEL_SIZE;
        goto out;
    }
    self->rng = gsl_rng_alloc(gsl_rng_default);
    if (self->rng == NULL) {
        ret = ERR_ALLOC_FAILED;
        goto out;
    }
    r = 0;
    for (j = 0; j < self->num_event_classes; j++) {
        if (self->event_classes[j].r > self->torus_diameter / 4) {
            ret = ERR_EVENT_TOO_LARGE;
            goto out;
        }
        if (self->event_classes[j].r > r) {
            r = self->event_classes[j].r;
        }
    }
    gsl_rng_set(self->rng, self->random_seed);
    self->P = xmalloc(gsl_pow_2(self->N) * sizeof(avl_tree_t));
    /* set up buffers */
    self->max_disc_pixels = gsl_pow_2(4 * ((int)(r / self->pixel_size) + 1));
    self->max_num_children = self->max_population_size;
    self->pixel_buffer = xmalloc(self->max_disc_pixels * sizeof(unsigned int));
    self->probability_buffer = xmalloc(self->max_occupancy * sizeof(double));
    self->intersected_buffer = xmalloc((self->max_occupancy + 1) 
            * sizeof(void *));
    self->child_buffer = xmalloc(self->max_num_children * sizeof(void *));
    self->parent_buffer = xmalloc((self->num_parents + 1) * sizeof(void *));
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
    assert(pixel < self->N * self->N);
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
    double r = self->event_classes[0].r;
    sim_get_disc_pixels(self, ind->location, r); 
    for (j = 1; j <= self->pixel_buffer[0]; j++) {
        pixel = self->pixel_buffer[j];
        ret = sim_add_individual_to_pixel(self, pixel, ind);
        if (ret != 0) {
            goto out;
        }
    }
    self->population_size++;
out:
    return ret;
}


static int 
sim_remove_individual_from_pixel(sim_t *self, unsigned int pixel, 
        individual_t *ind) 
{
    int ret = 0;
    avl_node_t *node;
    uintptr_t id = (uintptr_t) ind;
    uint64_t search = (uint64_t) id;
    unsigned int h;
    /* remove the id from this pixel */
    node = avl_search(&self->P[pixel], &search);
    assert(node != NULL);
    avl_unlink_node(&self->P[pixel], node);
    sim_free_avl_set_node(self, node);
    h = avl_count(&self->P[pixel]);
    ret = sim_add_pixel_to_occupancy(self, h, pixel);
    if (ret != 0) {
        goto out;
    }
    ret = sim_remove_pixel_from_occupancy(self, h + 1, pixel);
out:
    return ret;
}

static int 
sim_remove_individual(sim_t *self, individual_t *ind)
{
    int ret = 0;
    unsigned int j, pixel;
    double r = self->event_classes[0].r;
    sim_get_disc_pixels(self, ind->location, r); 
    for (j = 1; j <= self->pixel_buffer[0]; j++) {
        pixel = self->pixel_buffer[j];
        ret = sim_remove_individual_from_pixel(self, pixel, ind);
        if (ret != 0) {
            goto out;
        }
    }
    sim_free_individual(self, ind);
    self->population_size--;
out:
    return ret;

}
/*
 * Sets up the simulation so that we can enter the main loop.
 */
int
sim_initialise(sim_t *self)
{
    int ret = 0;
    unsigned int j, k, l, pixel;
    double r = self->event_classes[0].r;
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
        sim_get_disc_pixels(self, x, r); 
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
    self->population_size = self->sample_size;
    self->coalesced_loci[0] = 0;
    self->time = 0.0;
    self->num_reproduction_events = 0;
    self->num_non_reproduction_events = 0;
    self->ancestral_material = self->sample_size * self->num_loci;
    if (self->sample_size == 1) {
        /* when the sample size is 1, we simulate the history of the 
         * sample indefinitely
         */
        self->ancestral_material++;
    }
out:
    return ret;
}

static int
sim_get_large_event_children(sim_t *self, double *z, double r, double u)
{
    int ret = 0;
    uint64_t id;
    uintptr_t int_ptr;
    unsigned int pixel, j;
    individual_t *ind;
    avl_tree_t pop;
    avl_node_t *node, *pop_node;
    
    avl_init_tree(&pop, avl_set_compare, NULL); 
    sim_get_disc_pixels(self, z, r); 
    for (j = 1; j <= self->pixel_buffer[0]; j++) {
        pixel = self->pixel_buffer[j];
        for (node = self->P[pixel].head; node != NULL; node = node->next) {
            id = *((uint64_t *) node->item);
            int_ptr = (uintptr_t) id;
            ind = (individual_t *) int_ptr;
            if (torus_squared_distance(z, ind->location, self->torus_diameter) 
                    < gsl_pow_2(r)) {
                pop_node = sim_alloc_avl_set_node(self, id); 
                if (pop_node == NULL) {
                    ret = ERR_OUT_OF_AVL_SET_NODES;
                    goto out;
                }
                if (avl_insert_node(&pop, pop_node) == NULL) {
                    sim_free_avl_set_node(self, pop_node);
                }
            }
        }
    }
    /* Now go through the set of potential children and pick out those who 
     * were born.
     */
    self->num_children = 0;
    for (pop_node = pop.head; pop_node != NULL; pop_node = pop_node->next) {
        id = *((uint64_t *) pop_node->item);
        int_ptr = (uintptr_t) id;
        ind = (individual_t *) int_ptr;
        if (gsl_rng_uniform(self->rng) < u) {
            self->child_buffer[self->num_children] = ind;
            self->num_children++;
            if (self->num_children >= self->max_num_children) {
                ret = ERR_OUT_OF_INDIVIDUALS;
                goto out;
            }
        }
        sim_free_avl_set_node(self, pop_node);
    }
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
sim_generate_parents(sim_t *self)
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
    /* if we're simulating the pedigree, skip updating the ancestry */
    if (self->simulate_pedigree == 0) {
        for (j = 0; j < self->num_children; j++) {
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
    }
out:
    return ret;
}


/* 
 * Completes an event centred on z with radius r, assuming that the 
 * set of children has been collected into the child_buffer.
 */
static int
sim_complete_event(sim_t *self, double *z, double r)
{
    int ret = 0;
    unsigned int j, k;
    individual_t *ind;
    ret = sim_generate_parents(self);
    if (ret != 0) {
        goto out;
    }
    for (j = 0; j < self->num_children; j++) {
        ret = sim_remove_individual(self, self->child_buffer[j]);
        if (ret != 0) {
            goto out;
        }
    }
    /* add in the parents, if they have any ancestral material */
    for (j = 0; j < self->num_parents; j++) {
        ind = self->parent_buffer[j];
        /* We always add in both parents when simulating the pedigree */
        if (avl_count(&ind->ancestry) == 0 && self->simulate_pedigree == 0) {
            sim_free_individual(self, ind);
        } else {
            random_point_torus_disc(ind->location, z, r, 
                    self->torus_diameter, self->rng);
            if (self->dimension == 1) {
                ind->location[1] = 0.0;
            }
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
out:
    return ret;
}




/* 
 * Simulates a non-reproduction event.
 */
static int 
sim_non_reproduction_event(sim_t *self)
{
    int ret = 0;
    unsigned int j;
    double r, u;
    double z[] = {0.0, 0.0};
    double L = self->torus_diameter;
    double non_repr_rate = 0.0;
    double *p = xmalloc(self->num_event_classes * sizeof(double));
    
    for (j = 1; j < self->num_event_classes; j++) {
        non_repr_rate += self->event_classes[j].rate;
        p[j - 1] = self->event_classes[j].rate;
    }
    for (j = 0; j < self->num_event_classes - 1; j++) {
        p[j] /= non_repr_rate;
    }
    j = 1 + probability_list_select(p, self->num_event_classes - 1, 
            gsl_rng_uniform(self->rng));
    r = self->event_classes[j].r;
    u = self->event_classes[j].u;
    z[0] = L * gsl_rng_uniform(self->rng);
    z[1] = L * gsl_rng_uniform(self->rng);
    if (self->dimension == 1) {
        z[1] = 0.0;
    }
    ret = sim_get_large_event_children(self, z, r, u);
    if (ret != 0) {
        goto out;
    }
    ret = sim_complete_event(self, z, r);
    if (ret != 0) {
        goto out;
    }
    self->num_non_reproduction_events++;
out:
    free(p);
    return ret;
}


/*
 * Simulate the coalescent proces for at most the specified number 
 * of events. The simulation may terminate before this 
 * number of events occurs because the global time exceeds sim->max_time
 * or complete coalescence occurs. If the simulation finishes because 
 * of coalescence or time exceeded the stated maximum, 0 is returned.
 * Otherwise 1 is returned, indictating that more simulation is 
 * required. A negative return value indicates an error condition.
 */
int 
sim_simulate(sim_t *self, uint64_t max_events)
{
    int ret = 0;
    unsigned int j, k, pixel, max_occupancy, occupancy, pixel_count, v[2];
    double r = self->event_classes[0].r;
    double lambda = self->event_classes[0].rate;
    double non_repr_rate, total_rate;
    double Lambda_const = lambda * gsl_pow_2(self->pixel_size / self->torus_diameter);
    double z[] = {0.0, 0.0};
    double *p = self->probability_buffer;
    individual_t **S = self->intersected_buffer;
    unsigned int S_size;
    uint64_t set_item;
    uintptr_t int_ptr;
    double Lambda, jump_proba, *x;
    uint64_t events = 0;
    avl_node_t *node;
    set_map_value_t *smv, smv_search;
    individual_t *ind;

    if (self->dimension == 1) {
        /* In 1D we put all points on the line y=0. This means that we intersect
         * with exactly 2 pixels instead of 1, so we must divide the real constant 
         * by 2 to compensate for this.
         */
        assert(self->pixel_size == 2.0);
        Lambda_const = lambda * self->pixel_size / self->torus_diameter;
        Lambda_const /= 2;
    }
    non_repr_rate = 0.0;
    for (j = 1; j < self->num_event_classes; j++) {
        non_repr_rate += self->event_classes[j].rate;
    }
    /* we use Lambda = -1 to indicate that the state of the sample has 
     * changed and so we need to recalculate the rate of events */
    Lambda = -1.0;
    max_occupancy = 0;
    while (self->ancestral_material > self->num_loci && self->time < self->max_time
            && events < max_events) {
        /* first calculate the rate of (potential) reproduction events */
        if (Lambda == -1.0) {
            Lambda = 0.0;
            max_occupancy = 0;
            memset(p, 0, self->max_occupancy * sizeof(double)); 
            for (node = self->Q.head; node != NULL; node = node->next) {
                smv = (set_map_value_t *) node->item;
                occupancy = smv->key;
                pixel_count = avl_count(&smv->value);
                if (occupancy >= self->max_occupancy) {
                    ret = ERR_MAX_OCCUPANCY_EXCEEDED; 
                    goto out;
                }
                p[occupancy - 1] = pixel_count * sim_ubar(self, occupancy);
                Lambda += p[occupancy - 1]; 
                if (occupancy > max_occupancy) {
                    max_occupancy = occupancy;
                }
            }
            for (j = 1; j <= max_occupancy; j++) {
                p[j - 1] /= Lambda;
            }
        }
        /* Now determine the type of event this is */ 
        total_rate = Lambda * Lambda_const + non_repr_rate; 
        self->time += gsl_ran_exponential(self->rng, 1.0 / total_rate); 
        if (gsl_rng_uniform(self->rng) < non_repr_rate / total_rate) {
            ret = sim_non_reproduction_event(self); 
            if (ret != 0) {
                goto out;
            }
            events++;
            Lambda = -1.0;
        } else {
            /* this is a reproduction event */
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
            if (self->dimension == 1) {
                z[1] = 0.0;
            }
            assert(0.0 <= z[0] && self->torus_diameter > z[0]);
            assert(0.0 <= z[1] && self->torus_diameter > z[1]);
            S_size = 0;
            for (node = self->P[pixel].head; node != NULL; node = node->next) {
                set_item =  *((uint64_t *) node->item);
                int_ptr = (uintptr_t) set_item;
                ind = (individual_t *) int_ptr; 
                x = ind->location; 
                if (torus_squared_distance(z, x, self->torus_diameter) < gsl_pow_2(r)) {
                    S[S_size] = ind;
                    S_size++;
                }
            }
            jump_proba = sim_ubar(self, S_size) / sim_ubar(self, occupancy); 
            assert(jump_proba <= 1.0);
            if (gsl_rng_uniform(self->rng) < jump_proba) {
                events++;
                self->num_reproduction_events++;
                self->num_children = gsl_ran_discrete(self->rng, 
                        self->beta_distributions[S_size]);
                gsl_ran_choose(self->rng, self->child_buffer, self->num_children, 
                        S, S_size, sizeof(individual_t *));
                ret = sim_complete_event(self, z, r);
                if (ret != 0) {
                    goto out;
                }
                Lambda = -1.0;
            } 
        } 
    }
    ret = 1;
    if (self->ancestral_material == self->num_loci 
            || self->time >= self->max_time) {
        ret = 0;
    }
out:
    return ret;
}



int
sim_get_population(sim_t *self, avl_tree_t *pop)
{
    int ret = 0;
    uint64_t id;
    unsigned int pixel;
    avl_node_t *node, *pop_node;
    avl_init_tree(pop, avl_set_compare, NULL); 
    node = NULL; 
    for (pixel = 0; pixel < gsl_pow_2(self->N); pixel++) {
        for (node = self->P[pixel].head; node != NULL; node = node->next) {
            id = *((uint64_t *) node->item);
            pop_node = sim_alloc_avl_set_node(self, id); 
            if (pop_node == NULL) {
                ret = ERR_OUT_OF_AVL_SET_NODES;
                goto out;
            }
            if (avl_insert_node(pop, pop_node) == NULL) {
                sim_free_avl_set_node(self, pop_node);
            }
        }
    }
    assert(avl_count(pop) == self->population_size);
out: 
    return ret;
}

/* 
 * Frees the nodes in the specified set back to the simulation.
 */
void 
sim_free_population(sim_t *self, avl_tree_t *pop)
{
    avl_node_t *node;
    for (node = pop->head; node != NULL; node = node->next) {
        sim_free_avl_set_node(self, node);
    }
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

int
sim_print_state(sim_t *self, int detail)
{
    int ret = 0;
    unsigned int j, k, l, h, v[2];
    avl_tree_t chi;
    double *x;
    uint64_t id, pixel;
    uintptr_t int_ptr;
    set_map_value_t *set_map_value, search;
    avl_node_t *node, *chi_node;
    avl_init_tree(&chi, avl_set_compare, NULL); 
    individual_t *ind;
    node = NULL; 
    ret = sim_get_population(self, &chi);
    if (ret != 0) {
        goto out;
    }
    printf("reproduction events:    \t%lu\n", 
            (unsigned long) self->num_reproduction_events);
    printf("non reproduction events:\t%lu\n", 
            (unsigned long) self->num_non_reproduction_events);
    printf("chi = \n");
    for (chi_node = chi.head; chi_node != NULL; chi_node = chi_node->next) {
        id = *((uint64_t *) chi_node->item);
        int_ptr = (uintptr_t) id;
        ind = (individual_t *) int_ptr;
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
                printf(" %ld, ", (unsigned long) pixel);
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
    sim_free_population(self, &chi);
out:
    return ret;
}


