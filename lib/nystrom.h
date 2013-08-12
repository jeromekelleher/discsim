/*
 * Calculate the probability of identity in state of the disc model 
 * using the Nystrom method.
 *
 * author: Jerome Kelleher <jerome.kelleher@ed.ac.uk>
 */

#include <gsl/gsl_integration.h>

#ifdef __cplusplus
extern "C"{
#endif 

typedef struct nystrom_t_t {
    /* parameters */
    unsigned int nu;
    unsigned int n;
    double r;
    double u;
    double mu;
    double L;
    /* Nystrom solution vectors */
    double *x;
    double *w;
    double *f;
    /* GSL integration parameters */
    size_t integration_workspace_size;
    gsl_integration_workspace *integration_workspace;
    int integration_rule;
    double integration_abserr;
    double integration_relerr;
} nystrom_t;

/*
 * Allocate a new nystrom solver for the probability of identity in state.
 * Events with the specified properties occur at rate 1 on a torus of diameter
 * L. See 
 * http://www.gnu.org/software/gsl/manual/html_node/QAG-adaptive-integration.html
 * for details of the parameters used for numerical integration.
 * 
 * :param r: the radius of events
 * :param u: the intensity of events
 * :param mu: the mutation rate
 * :param nu: the number of parents in an event
 * :param L: the diameter of the torus (this just sets the rate)
 * :param max_x: the maximum x that we solve for
 * :param n: the number of quadrature points used to solve the Fredholm 
 *      equation.
 * :param integration_workspace_size: the size of the GSL integration workspace
 *      used for numerical integration
 * :param integration_abserr: the desired absolute error for numerical 
 *     integration.
 * :param integration_relerr: the desired relative error for numerical 
 *     integration. See the GSL documentation for a description of these
 *     error bounds.
 */
nystrom_t * nystrom_alloc(double r, double u, double mu, int nu, double L,
        double max_x, unsigned int n, unsigned int integration_workspace_size,
        double integration_abserr, double integration_relerr);

/* 
 * Frees the memory associated with the specified nystrom solver.
 */
void nystrom_free(nystrom_t *self);

/* 
 * Solves the integral equation and make the object ready for interpolation.
 */
void nystrom_solve(nystrom_t *self);

/*
 * Returns the value of F(x) using Nystrom interpolation. 
 */
double nystrom_interpolate(nystrom_t *self, double x);

/* debug method */
void nystrom_print_state(nystrom_t *self);

#ifdef __cplusplus
}
#endif
