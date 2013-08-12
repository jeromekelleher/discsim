/*
 * Calculate the probability of identity in state of the disc model 
 * using the Nystrom method.
 *
 * author: Jerome Kelleher <jerome.kelleher@ed.ac.uk>
 */

#include <math.h>
#include <assert.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>

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
    double r = self->r;
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
    double r = self->r;
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
    double phi = 2 * self->mu + 2 * self->u * M_PI * gsl_pow_2(self->r) 
               - gsl_pow_2(self->u) * A2(x, self->r);
    return phi;
}

static double 
nystrom_K(nystrom_t *self, double x, double R)
{
    double r = self->r;
    double u = self->u;
    int nu = self->nu;
    double K = u * u * disc_distance_pdf(R / r) * (1.0 - 1.0 / nu) * A2(x, r) / r 
            + 4 * u * R * (intQ1(self, x, R) - u * intQ2(self, x, R))
            / (M_PI * gsl_pow_2(r));
    return K;
}

static double 
nystrom_g(nystrom_t *self, double x)
{
    return gsl_pow_2(self->u) * A2(x, self->r) / self->nu;
}    

/* 
 * public interface 
 */
void
nystrom_print_state(nystrom_t *self) 
{
    unsigned int j;
    printf("r = %f\n", self->r);
    printf("u = %f\n", self->u);
    printf("mu = %f\n", self->mu);
    printf("L = %f\n", self->L);
    printf("nu = %d\n", self->nu);
    printf("n = %d\n", self->n);
    printf("x = ");
    for (j = 0; j < self->n; j++) {
        printf("%f, ", self->x[j]);
    }
    printf("\nw = ");
    for (j = 0; j < self->n; j++) {
        printf("%f, ", self->w[j]);
    }
    printf("\n");
}


nystrom_t *
nystrom_alloc(double r, double u, double mu, int nu, double L, double max_x,
        unsigned int n, unsigned int integration_workspace_size,
        double integration_abserr, double integration_relerr)
{
    nystrom_t *self = malloc(sizeof(nystrom_t));
    unsigned int i;
    double a = 0.0;
    double b = max_x;
    double A = (b - a) / 2.0;
    double B = (b + a) / 2.0;
    gsl_integration_glfixed_table *t = gsl_integration_glfixed_table_alloc(n);
    self->r = r;
    self->u = u;
    self->mu = mu * L * L;
    self->nu = nu;
    self->L = L;
    self->n = n;
    self->w = malloc(n * sizeof(double));
    self->x = malloc(n * sizeof(double));
    self->f = malloc(n * sizeof(double));
    for (i = 0; i < n / 2; i++) {
        self->w[i] = t->w[n / 2 - i - 1] * A;
        self->w[n - i - 1] = self->w[i]; 
        self->x[i] = B - A * t->x[n / 2 - i - 1];
        self->x[n - i - 1] = B + A * t->x[n / 2 - i - 1];
    }
    gsl_integration_glfixed_table_free(t);
    /* Set up the numerical integration parameters */
    self->integration_rule = GSL_INTEG_GAUSS31; 
    self->integration_workspace_size = integration_workspace_size;
    self->integration_abserr = integration_abserr;
    self->integration_relerr = integration_relerr; 
    self->integration_workspace = gsl_integration_workspace_alloc(
            self->integration_workspace_size);
    return self;
}

void
nystrom_free(nystrom_t *self)
{
    free(self->w);
    free(self->x);
    free(self->f);
    gsl_integration_workspace_free(self->integration_workspace);
    free(self);
}

void
nystrom_solve(nystrom_t *self)
{
    unsigned int i, j;
    size_t n = self->n;
    double *x = self->x;
    double *w = self->w;
    double K, phi;
    int s;
    gsl_permutation *p = gsl_permutation_alloc(n);
    gsl_matrix *kernel = gsl_matrix_alloc(n, n);
    gsl_vector *gv = gsl_vector_alloc(n);
    gsl_vector_view fvv = gsl_vector_view_array(self->f, n);
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
    gsl_permutation_free(p);
    gsl_matrix_free(kernel);
    gsl_vector_free(gv);
}


double 
nystrom_interpolate(nystrom_t *self, double x)
{
    unsigned int i;
    double K;
    double sum = 0.0;
    double phi = nystrom_phi(self, x);
    for (i = 0; i < self->n; i++) {
        K = nystrom_K(self, x, self->x[i]) / phi;
        sum += self->w[i] * self->f[i] * K; 
    }
    return nystrom_g(self, x) / phi + sum;
}



