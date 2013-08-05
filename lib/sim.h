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

#include <stdint.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


#define ERR_ALLOC_FAILED 1
#define ERR_BAD_PIXEL_SIZE 2 
#define ERR_OUT_OF_INT_MAP_NODES 3
#define ERR_OUT_OF_SET_MAP_NODES 4 
#define ERR_OUT_OF_AVL_SET_NODES 5 
#define ERR_OUT_OF_INDIVIDUALS 6 
#define ERR_AVL_OP_FAILED 7

#include "avl.h"

typedef struct {
    double location[2];
    avl_tree_t ancestry;
} individual_t;

typedef struct {
    unsigned int key;
    avl_tree_t value;
} set_map_value_t;

typedef struct {
    unsigned int key;
    unsigned int value;
} int_map_value_t;

typedef struct {
    double rate;
    double r;
    double u;
} event_class_t;

typedef struct {
    double *sample;
    unsigned int sample_size;
    unsigned int num_loci;
    unsigned int num_parents;
    unsigned int max_population_size;
    unsigned int max_occupancy;
    unsigned long random_seed;
    double torus_diameter;
    double pixel_size;
    double recombination_probability;
    /* Events */ 
    unsigned int num_event_classes;
    event_class_t *event_classes;
    /* algorithm state */
    double beta_threshold;
    unsigned int max_disc_pixels;
    unsigned int N;
    gsl_rng *rng;
    avl_tree_t *P;
    avl_tree_t Q;
    double time;
    uint64_t successful_events;
    unsigned int population_size;
    unsigned int ancestral_material;
    double *ubar;
    gsl_ran_discrete_t **beta_distributions;
    /* set memory */ 
    uint64_t *avl_set_value_mem;
    avl_node_t *avl_set_node_mem;
    avl_node_t **avl_set_node_heap;
    int avl_set_node_heap_top;
    /* individual memory */ 
    individual_t *individual_mem;
    individual_t **individual_heap;
    int individual_heap_top;
    /* set map memory */ 
    set_map_value_t *avl_set_map_value_mem;
    avl_node_t *avl_set_map_node_mem;
    avl_node_t **avl_set_map_node_heap;
    int avl_set_map_node_heap_top;
    /* int map memory */ 
    int_map_value_t *avl_int_map_value_mem;
    avl_node_t *avl_int_map_node_mem;
    avl_node_t **avl_int_map_node_heap;
    int avl_int_map_node_heap_top;
    /* buffers */
    unsigned int *pixel_buffer;
    double *probability_buffer;
    individual_t **intersected_buffer;
    individual_t **child_buffer;
    individual_t **parent_buffer;
    /* ancestry */
    int **pi;
    double **tau;
    int *eta;
    /* memory for the large ancestry blocks */
    int *pi_mem;
    double *tau_mem;
    /* Map to keep track of coalescences */
    unsigned int *coalescence_map;
    unsigned int *coalesced_loci;
} sim_t;

void sim_free(sim_t *self);
void sim_set_max_occupancy(sim_t *self, double headroom);
int sim_alloc(sim_t *self);
int sim_initialise(sim_t *self);
int sim_simulate(sim_t *self, unsigned int max_jumps, double max_time);
void sim_print_parameters(sim_t *self);
void sim_print_state(sim_t *self, int detail);
const char * sim_error_message(int err);

