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
    unsigned int n;
    unsigned int m;
    unsigned int nu;
    double L;
    double s;
    double r;
    double u;
    double lambda;
    double *X;
    double rho;
    unsigned long random_seed;
    unsigned int max_sample_size;
    unsigned int max_occupancy;
    double max_time;
    double beta_threshold;
    /* algorithm state */
    unsigned int max_disc_pixels;
    unsigned int N;
    gsl_rng *rng;
    avl_tree_t *P;
    avl_tree_t Q;
    double t;
    unsigned int sample_size;
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


