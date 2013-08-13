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

#include <gsl/gsl_integration.h>

typedef struct {
    /* parameters */
    unsigned int num_quadrature_points;
    unsigned int num_parents;
    unsigned int num_event_classes;
    event_class_t *event_classes;
    double mutation_rate;
    double torus_diameter;
    double max_x;
    /* GSL integration parameters */
    unsigned int integration_workspace_size;
    double integration_abserr;
    double integration_relerr;
    int integration_rule;
    /* Internal state */
    double *x;
    double *w;
    double *f;
    gsl_integration_workspace *integration_workspace;
} nystrom_t;

int nystrom_alloc(nystrom_t *self);
/* 
 * Frees the memory associated with the specified nystrom solver.
 */
void nystrom_free(nystrom_t *self);

/* 
 * Solves the integral equation and make the object ready for interpolation.
 */
int nystrom_solve(nystrom_t *self);

/*
 * Returns the value of F(x) using Nystrom interpolation. 
 */
double nystrom_interpolate(nystrom_t *self, double x);

/* debug method */
void nystrom_print_state(nystrom_t *self);

