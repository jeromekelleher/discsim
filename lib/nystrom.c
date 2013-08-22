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

#include <math.h>
#include <assert.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>

#include "util.h"
#include "nystrom.h"

typedef struct {
    nystrom_t *nystrom;
    double x;
} K_integral_params;

typedef struct {
    double x;
    double R;
    double r;
} Q_integral_params_t;

/* 
 * Returns the probability density of the distance x between two points 
 * selected uniformly at random within the unit disc.
 */
static double 
disc_distance_pdf(double x) 
{
    double f = 0.0;
    if (0 <= x && x < 2.0) {
        f = x * (4 * acos(x / 2) - x * sqrt(4 - gsl_pow_2(x))) / M_PI;
    }
    return f;
}

/*
 * Returns the area of the intersection of two circles with radius r 
 * with centres distance d apart.
 */
static double 
A2(double d, double r)
{
    double ret = 0.0;
    if (d < 2 * r) {
        ret = 2 * r * r * acos(d / (2 * r)) - d * sqrt(4 * r * r - d * d) / 2;
    }
    return ret;
}

static inline double
euclidean_distance(double x[2], double y[2])
{
    return sqrt(gsl_pow_2(x[0] - y[0]) + gsl_pow_2(x[1] - y[1]));
}

static double 
cta(double d, double z[2], double r)
{
    unsigned int j;
    double h = sqrt(gsl_pow_2(r) - gsl_pow_2(d) / 4);
    double c[3];
    double p1[2];
    double p2[2];
    double p3[2];
    double A;
    p1[0] = (-d + 2*(z[0] + (z[1]*sqrt(16*pow(r,2) - pow(d + 2*z[0],2) 
            - 4*pow(z[1],2)))/sqrt(pow(d + 2*z[0],2) + 4*pow(z[1],2))))/4;
    p1[1] = z[1]/2. - ((d + 2*z[0])*sqrt(16*pow(r,2) - pow(d + 2*z[0],2) 
            - 4*pow(z[1],2)))/(4*sqrt(pow(d + 2*z[0],2) + 4*pow(z[1],2)));
    p2[0] = (d + 2*z[0] - (2*z[1]*sqrt(16*pow(r,2) - pow(d - 2*z[0],2) 
            - 4*pow(z[1],2)))/sqrt(pow(d - 2*z[0],2) + 4*pow(z[1],2)))/4.;
    p2[1] = z[1]/2. - ((d - 2*z[0])*sqrt(16*pow(r,2) - pow(d - 2*z[0],2) 
            - 4*pow(z[1],2)))/(4.*sqrt(pow(d - 2*z[0],2) + 4*pow(z[1],2)));
    p3[0] = 0.0;
    p3[1] = h;
    c[0] = euclidean_distance(p1, p2);
    c[1] = euclidean_distance(p2, p3);
    c[2] = euclidean_distance(p1, p3);
    A = sqrt((c[0] + c[1] + c[2]) * (-c[0] + c[1] + c[2]) * (c[0] - c[1] + c[2])
            * (c[0] + c[1] - c[2])) / 4;
    for (j = 0; j < 3; j++) {
        A += gsl_pow_2(r) * asin(c[j] / (2 * r)) - c[j] 
                * sqrt(4 * gsl_pow_2(r) - gsl_pow_2(c[j])) / 4;
    }
    return A;
}

/* 
 * Returns the area of the intersection of three circles, centered at 
 * (-d/2,0), (d/2,0) and z. 
 */
static double 
A3(double d, double y[2], double r)
{
    double ret = 0.0;
    double h = sqrt(gsl_pow_2(r) - gsl_pow_2(d) / 4);
    double r2 = gsl_pow_2(r);
    double z[2];
    /* map z to the upper quarter plane */
    z[0] = fabs(y[0]);
    z[1] = fabs(y[1]);
    if (gsl_pow_2(z[0]) + gsl_pow_2(z[1] + h) <= r2) {
        /* 
         * if z is within a circle of radius r centered at (0, -h) the 
         * third disc entirely overlaps.
         */
        ret = A2(d, r);
    } else if (gsl_pow_2(z[0]) + gsl_pow_2(z[1] - h) <= r2) {
        /*
         * Otherwise, if z is within a circle of radius r centered at (0, h)
         * we have a circular triangle.
         */
        ret = cta(d, z, r);
    } else if (gsl_pow_2(z[0] + d / 2) + gsl_pow_2(z[1]) <= 4 * r2) {
        /*
         * Otherwise, if z is within a circle or radius 2r centered at 
         * (-d/2, 0) we have a lens overlap between the circle at (-d/2,0)
         * and z.
         */
        ret = A2(sqrt(gsl_pow_2(z[0] + d / 2) + gsl_pow_2(z[1])), r);
    }
    /* otherwise, the intersection is empty */
    return ret;
}


static double 
intQ1_f(double theta, void *params) {
    Q_integral_params_t *p = (Q_integral_params_t *) params;
    double x = p->x;
    double R = p->R;
    double r = p->r;
    /* we must protect against -0.0 here - sigh */
    double d2 = fabs(gsl_pow_2(x) - 2 * x * R * cos(theta) + gsl_pow_2(R));
    double ret = A2(sqrt(d2), r);
    return ret;
}

static double 
intQ2_f(double theta, void *params) {
    Q_integral_params_t *p = (Q_integral_params_t *) params;
    double x = p->x;
    double R = p->R;
    double r = p->r;
    double z[2];
    z[0] = -0.5 * x + R * cos(theta);
    z[1] = R * sin(theta);
    return A3(x, z, r);
}

static double 
intQ2(nystrom_t *self, double x, double R)
{
    double ret = 0.0;
    double v, error, alpha, beta, gamma, d;
    double r = self->event_classes[0].r;
    double r2 = gsl_pow_2(r);
    double x2 = gsl_pow_2(x);
    double R2 = gsl_pow_2(R);
    double s = sqrt((R2 - 4 * r2) * (x2 - 4 * r2));
    double nu, kappa;
    Q_integral_params_t p;
    gsl_function f;
    p.x = x;
    p.R = R;
    p.r = r;
    f.function = &intQ2_f;
    f.params = &p;
    if (x == 0.0) {
        ret = M_PI * A2(R, r);
    } else if (R == 0.0) {
        ret = M_PI * A2(x, r);
    } else if (x < 2 * r && R < 2 * r) { 
        nu =  (R * x + s) / (4 * r2);
        alpha = -1.0 < nu && nu < 1.0 ? acos(nu): 0.0;
        nu =  (R * x - s) / (4 * r2);
        assert(-1.0 < nu && nu < 1.0);
        beta = acos((R * x - s) / (4 * r2));
        /* Take the initial constant section */
        d = R < x? x : R;
        ret = alpha * A2(d, r);
        gsl_integration_qag(&f, alpha, beta, self->integration_abserr, 
                self->integration_relerr, self->integration_workspace_size, 
                self->integration_rule, self->integration_workspace, 
                &v, &error); 
        ret += v;
        if (beta > M_PI / 2) {
            gamma = M_PI;
            kappa = (R2 - 4 * r2 + x2) / (2 * R * x);
            if (-1.0 < kappa && 1.0 > kappa) {
                gamma = acos(kappa);
            }
            f.function = &intQ1_f;
            gsl_integration_qag(&f, beta, gamma, self->integration_abserr, 
                    self->integration_relerr, self->integration_workspace_size,
                    self->integration_rule, self->integration_workspace, 
                    &v, &error); 
            ret += v;
        }
    }
    return ret; 
}

static double 
intQ1(nystrom_t *self, double x, double R)
{
    double ret = 0.0;
    double error, alpha, nu;
    double r = self->event_classes[0].r;
    Q_integral_params_t p;
    gsl_function f;
    p.x = x;
    p.R = R;
    p.r = r;
    f.function = &intQ1_f;
    f.params = &p;
    if (x == 0.0) {
        ret = M_PI * A2(R, r);
    } else if (R == 0.0) {
        ret = M_PI * A2(x, r);
    } else {
        alpha = M_PI;
        nu = (gsl_pow_2(x) + gsl_pow_2(R) - 4 * gsl_pow_2(r)) / (2 * x * R);
        if (-1.0 < nu && nu < 1.0) { 
            alpha = acos(nu);
        }
        gsl_integration_qag(&f, 0.0, alpha, self->integration_abserr, 
                self->integration_relerr, self->integration_workspace_size, 
                self->integration_rule, self->integration_workspace, 
                &ret, &error); 
    }
    return ret; 
}

static double 
nystrom_phi(nystrom_t *self, double x)
{
    /* translate rates into torus rates */
    double lambda = self->event_classes[0].rate;
    double Lambda = lambda / gsl_pow_2(self->torus_diameter);
    double mu = self->mutation_rate / Lambda; 
    double u = self->event_classes[0].u;
    double r = self->event_classes[0].r;
    double phi = 2 * mu + 2 * u * M_PI * gsl_pow_2(r) 
            - gsl_pow_2(u) * A2(x, r);
    return phi;
}

static double 
nystrom_K(nystrom_t *self, double x, double R)
{
    double r = self->event_classes[0].r;
    double u = self->event_classes[0].u;
    int nu = self->num_parents;
    double K = u * u * disc_distance_pdf(R / r) * (1.0 - 1.0 / nu) * A2(x, r) / r 
            + 4 * u * R * (intQ1(self, x, R) - u * intQ2(self, x, R))
            / (M_PI * gsl_pow_2(r));
    return K;
}

static double 
nystrom_g(nystrom_t *self, double x)
{
    double r = self->event_classes[0].r;
    double u = self->event_classes[0].u;
    double nu = self->num_parents;
    return gsl_pow_2(u) * A2(x, r) / nu;
}    

/* 
 * public interface 
 */
void
nystrom_print_state(nystrom_t *self) 
{
    unsigned int j;
    printf("r = %f\n", self->event_classes[0].r);
    printf("u = %f\n", self->event_classes[0].u);
    printf("mu = %f\n", self->mutation_rate);
    printf("L = %f\n", self->torus_diameter);
    printf("nu = %d\n", self->num_parents);
    printf("n = %d\n", self->num_quadrature_points);
    printf("x = ");
    for (j = 0; j < self->num_quadrature_points; j++) {
        printf("%f, ", self->x[j]);
    }
    printf("\nw = ");
    for (j = 0; j < self->num_quadrature_points; j++) {
        printf("%f, ", self->w[j]);
    }
    printf("\n");
}


int
nystrom_alloc(nystrom_t *self)
{
    int ret = 0;
    unsigned int j;
    unsigned int n = self->num_quadrature_points;
    double a = 0.0;
    double b = self->max_x;
    double A = (b - a) / 2.0;
    double B = (b + a) / 2.0;
    gsl_integration_glfixed_table *t = gsl_integration_glfixed_table_alloc(n);
    /* this can easily be made a parameter */ 
    self->integration_rule = GSL_INTEG_GAUSS31; 
    self->integration_workspace = gsl_integration_workspace_alloc(
            self->integration_workspace_size);
    self->w = malloc(n * sizeof(double));
    self->x = malloc(n * sizeof(double));
    self->f = malloc(n * sizeof(double));
    if (t == NULL || self->x == NULL || self->w == NULL || self->f == NULL
            || self->integration_workspace == NULL) {
        ret = ERR_ALLOC_FAILED;
        goto out;
    }
    for (j = 0; j < n / 2; j++) {
        self->w[j] = t->w[n / 2 - j - 1] * A;
        self->w[n - j - 1] = self->w[j]; 
        self->x[j] = B - A * t->x[n / 2 - j - 1];
        self->x[n - j - 1] = B + A * t->x[n / 2 - j - 1];
    }
    for (j = 0; j < n; j++) {
        self->f[j] = 0.0;
    }   
out:
    return ret;
}

void
nystrom_free(nystrom_t *self)
{
    free(self->w);
    free(self->x);
    free(self->f);
    gsl_integration_workspace_free(self->integration_workspace);
}

int
nystrom_solve(nystrom_t *self)
{
    int ret = 0;
    unsigned int i, j;
    size_t n = self->num_quadrature_points;
    double *x = self->x;
    double *w = self->w;
    double K, phi;
    int s;
    gsl_vector_view fvv = gsl_vector_view_array(self->f, n);
    gsl_permutation *p = NULL;
    gsl_matrix *kernel = NULL; 
    gsl_vector *gv = NULL; 
   
    p = gsl_permutation_alloc(n);
    kernel = gsl_matrix_alloc(n, n);
    gv = gsl_vector_alloc(n);
    if (p == NULL || kernel == NULL || gv == NULL) {
        ret = -ERR_ALLOC_FAILED;
        goto out;
    }
    for (i = 0; i < n; i++) {
        phi = nystrom_phi(self, x[i]);
        for (j = 0; j < n; j++) {
            K = nystrom_K(self, x[i], x[j]) / phi;
            gsl_matrix_set(kernel, i, j, (double) (i == j) - w[j] * K); 
        }
        gsl_vector_set(gv, i, nystrom_g(self, x[i]) / phi);
    }
    gsl_linalg_LU_decomp(kernel, p, &s);
    gsl_linalg_LU_solve(kernel, p, gv, &(fvv.vector));
out:
    if (p != NULL) {
        gsl_permutation_free(p);
    }
    if (kernel != NULL) {
        gsl_matrix_free(kernel);
    }
    if (gv != NULL) {
        gsl_vector_free(gv);
    }
    return ret;
}


double 
nystrom_interpolate(nystrom_t *self, double x)
{
    unsigned int j;
    unsigned int n = self->num_quadrature_points;
    double K;
    double sum = 0.0;
    double phi = nystrom_phi(self, x);
    for (j = 0; j < n; j++) {
        K = nystrom_K(self, x, self->x[j]) / phi;
        sum += self->w[j] * self->f[j] * K; 
    }
    return nystrom_g(self, x) / phi + sum;
}



