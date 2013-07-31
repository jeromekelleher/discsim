#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <sys/times.h>
#include <unistd.h>
#include <stdint.h>
#include <avl.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_histogram.h>

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
    /* counters for analysis */
    unsigned long num_generated_events;
    unsigned long num_jumps;
    unsigned long max_jumps;
    unsigned long *total_E_size;
    unsigned long *total_H_size;
    unsigned long *total_S_size;
    unsigned long *total_C_size;
    unsigned long *total_G_size;
    gsl_histogram *hit_rate_distribution;
} sim_t;

typedef struct {
    unsigned int n;
    double L;
    double s;
    double r;
    double u;
    double lambda;
    double *X;
    double mu;
    double max_time;
    unsigned long random_seed;
    /* algorithm state */
    unsigned int N;
    avl_tree_t *P;
    avl_tree_t *Q;
    unsigned int max_occupancy;
    unsigned int max_disc_pixels;
    /* avl node memory */
    unsigned int *avl_intset_value_mem;
    avl_node_t *avl_intset_node_mem;
    avl_node_t **avl_intset_node_heap;
    int avl_intset_node_heap_top;
    /* output */
    int *pi;
    double *tau;
    double *chi;
    gsl_rng *rng;
    /* counters for analysis */
    unsigned long num_generated_events;
    unsigned long total_max_occupancy;
    unsigned long num_jumps;
    unsigned long max_jumps;
} alg_N_sim_t;

typedef struct {
    int *lambda;
    int *pi;
    int *tau;
    int *beta;
    int *alpha;
    void *__mem;
} sv_t;


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
 * Nearest Common Ancestor computation. An implementation of Schieber and Vishkin's 
 * nearest common ancestor algorithm from TAOCP 7.1.3 pg.34-35. Preprocesses the 
 * input Oriented Forest in O(n) time and processes queries for the nearest common 
 * ancestor between an arbitary pair of nodes in O(1) time.
 */

#define LAMBDA 0
/*
 * Computes the tables required to compute the nearest common ancestor using the 
 * Schieber-Vishkin algorithm for the specified set of vertices (numbered 1 to n, 
 * including the null vertex, LAMBDA) and the specified parent array, where v = p[u]
 * means that v is the parent of int u. 
 */  
void
sv_init(sv_t *self, int num_vertices, int *arcs) 
{
    int u, v, p, n, h, *ptr, *local_mem, not_done, *child, *parent, *sib;
    /* allocate the sv_tables and storage for it. */
    self->__mem = NULL;
    self->__mem = xcalloc((size_t) (5 * num_vertices), sizeof(int));
    ptr = self->__mem;
    self->lambda = ptr;
    ptr += num_vertices;
    self->pi = ptr; 
    ptr += num_vertices;
    self->tau = ptr; 
    ptr += num_vertices;
    self->alpha = ptr; 
    ptr += num_vertices;
    self->beta = ptr; 
    
    local_mem = xmalloc(3 * (size_t) num_vertices * sizeof(int));
    ptr = local_mem; 
    /* allocate the temporary local storage for navigating the tree. */
    child = ptr; 
    ptr += num_vertices;
    parent = ptr;
    ptr += num_vertices;
    sib = ptr; 
  
    /* set up the triply linked tree structure. */
    for (v = 0; v < num_vertices; v++) {
        child[v] = LAMBDA;
        sib[v] = LAMBDA;
        parent[v] = LAMBDA;
    }
    for (u = 1; u < num_vertices; u++) {
        v = arcs[u];
        sib[u] = child[v];
        child[v] = u;
        parent[u] = v;
    }

    p = child[LAMBDA];
    n = 0;
    self->lambda[0] = -1;
    
    while (p != LAMBDA) {
        not_done = 1;
        while (not_done) {
            n++;
            self->pi[p] = n;
            self->tau[n] = LAMBDA;
            self->lambda[n] = 1 + self->lambda[n >> 1];
            if (child[p] != LAMBDA) {
                p = child[p];
            } else {
                not_done = 0;
            }
        }    
        self->beta[p] = n;
        not_done = 1; 
        while (not_done) {
            self->tau[self->beta[p]] = parent[p];
            if (sib[p] != LAMBDA) {
                p = sib[p];
                not_done = 0;
            } else {
                p = parent[p];
                if (p != LAMBDA) {
                    h = self->lambda[n & -self->pi[p]];
                    self->beta[p] = ((n >> h) | 1) << h;
                } else {
                    not_done = 0;
                }
            }     
        }
    }

    /* Begin the second traversal */
    self->lambda[0] = self->lambda[n];
    self->pi[LAMBDA] = 0;
    self->beta[LAMBDA] = 0;
    self->alpha[LAMBDA] = 0;
    p = child[LAMBDA];
    
    while (p != LAMBDA) {
        not_done = 1;
        while (not_done) {
            self->alpha[p] = self->alpha[parent[p]] | (self->beta[p] & -self->beta[p]);
            if (child[p] != LAMBDA) {
                p = child[p];
            } else {
                not_done = 0;
            }
        }    
        not_done = 1;
        while (not_done) {
            if (sib[p] != LAMBDA) {
                p = sib[p];
                not_done = 0;
            } else {
                p = parent[p];
                not_done = p != LAMBDA;
            }
        }
    }    
    free(local_mem);
}

/*
 * Returns the nearest common ancestor of the nodes x and y.
 */
inline int 
sv_nca(sv_t *self, int x, int y)
{
    int z, xhat, yhat, k, h, j, l;

    if (self->beta[x] <= self->beta[y]) {
        h = self->lambda[self->beta[y] & -self->beta[x]];
    } else {
        h = self->lambda[self->beta[x] & -self->beta[y]];
    }    
    k = self->alpha[x] & self->alpha[y] & -(1 << h);
    h = self->lambda[k & -k];
    j = ((self->beta[x] >> h) | 1) << h;
    if (j == self->beta[x]) {
        xhat = x;
    } else {
        l = self->lambda[self->alpha[x] & ((1 << h) - 1)];
        xhat = self->tau[((self->beta[x] >> l) | 1) << l];
    }    
    if (j == self->beta[y]) {
        yhat = y;
    } else {
        l = self->lambda[self->alpha[y] & ((1 << h) - 1)];
        yhat = self->tau[((self->beta[y] >> l) | 1) << l];
    }    
    if (self->pi[xhat] <= self->pi[yhat]) {
        z = xhat;
    } else {
        z = yhat;
    }     
    return z; 
}

void
sv_free(sv_t *self) 
{
    free(self->__mem);
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



/* 
 * Runs the pairwise coalescent algorithm and returns the resulting 
 * coalcensnce time.
 */
double 
sim_algorithm_P(sim_t *self, double x)
{
    int j, k;
    double chi[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
    double z[2] = {0.0, 0.0};
    double t = 0.0;
    double d, Lambda, A1, A2, p; //d0, d1;
    double u = self->u;
    double r2 = gsl_pow_2(self->r);
    double L2 = gsl_pow_2(self->L);
    int not_done = 1;
    chi[0][0] = x;
    while (not_done) {
        d = sqrt(torus_squared_distance(chi[0], chi[1], self->L));
        A2 = disc_intersection_area(d, self->r);
        A1 = 2 * M_PI * r2 -  2 * A2;
        Lambda = self->lambda * (A1 * u + A2 * (1.0 - gsl_pow_2(1.0 - u))) / L2;
        t += gsl_ran_exponential(self->rng, 1.0 / Lambda);
        p = A1 / (A1 - A2 * (u - 2));
        j = gsl_rng_uniform(self->rng) < 0.5 ? 0 : 1;
        k = (j + 1)  % 2;
        if (gsl_rng_uniform(self->rng) < p) {
            /* Choose z from symmetric difference */
            assert(d > 0);
            do {
                random_point_torus_disc(z, chi[j], self->r, self->L, self->rng);
            } while (torus_squared_distance(z, chi[k], self->L) <= r2);
        } else {
            /* Choose z from intersection */
            assert(d < 2 * self->r);
            do {
                random_point_torus_disc(z, chi[j], self->r, self->L, self->rng);
            } while (torus_squared_distance(z, chi[k], self->L) > r2);
            if (gsl_rng_uniform(self->rng) <  u / (2 - u)) {
                not_done = 0;
            }
        }
        random_point_torus_disc(chi[j], z, self->r, self->L, self->rng);
        
    }
    return t;
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



/*=============================================================================
 *
 * Algorithm M implementation 
 *
 *=============================================================================
 */

void 
sim_write_occupancy_distribution(sim_t *self, char *filename) 
{
    unsigned int j, o;
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        fatal_error("IO error opening file");
    }
    for (j = 0; j < gsl_pow_2(self->N); j++) {
        o = avl_count(&self->P[j]);
        fprintf(f, "%d\t", o); 
        if ((j + 1) % self->N == 0) {
            fprintf(f, "\n");
        }       
    }
    fclose(f);
}

void 
sim_write_ancestry_distribution(sim_t *self, char *filename) 
{
    unsigned int j;
    unsigned int ancestral_material; 
    uint64_t set_value;
    individual_t *ind;
    FILE *f = fopen(filename, "w");
    avl_node_t *node;
    if (f == NULL) {
        fatal_error("IO error opening file");
    }
    for (j = 0; j < gsl_pow_2(self->N); j++) {
        ancestral_material = 0;
        for (node = self->P[j].head; node != NULL; node = node->next) {
            set_value = *((uint64_t *) node->item);
            ind = (individual_t *) set_value;
            ancestral_material += avl_count(&ind->ancestry);
        }
        fprintf(f, "%d\t", ancestral_material); 
        if ((j + 1) % self->N == 0) {
            fprintf(f, "\n");
        }       
    }
    fclose(f);
}

void 
sim_write_block_distribution(sim_t *self, char *filename) 
{
    unsigned int j;
    unsigned int last_locus, num_blocks;
    int_map_value_t *locus_mapping;
    uint64_t set_value;
    individual_t *ind;
    FILE *f = fopen(filename, "w");
    avl_node_t *node, *anode;
    if (f == NULL) {
        fatal_error("IO error opening file");
    }
    for (j = 0; j < gsl_pow_2(self->N); j++) {
        num_blocks = 0;
        for (node = self->P[j].head; node != NULL; node = node->next) {
            set_value = *((uint64_t *) node->item);
            ind = (individual_t *) set_value;
            anode = ind->ancestry.head;
            assert(anode != NULL);
            locus_mapping = (int_map_value_t *) anode->item;
            last_locus = locus_mapping->key;
            num_blocks++;
            while (anode->next != NULL) {
                anode = anode->next; 
                assert(anode != NULL);
                locus_mapping = (int_map_value_t *) anode->item;
                if (locus_mapping->key != last_locus + 1) {
                    num_blocks++;
                }
                last_locus = locus_mapping->key;
            }
        }
        fprintf(f, "%d\t", num_blocks); 
        if ((j + 1) % self->N == 0) {
            fprintf(f, "\n");
        }       
    }
    fclose(f);
}

void 
sim_write_distribution(sim_t *self, char *filename, unsigned long *matrix)
{
    unsigned int j;
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        fatal_error("IO error opening file");
    }
    for (j = 0; j < gsl_pow_2(self->N); j++) {
        fprintf(f, "%ld\t", matrix[j]); 
        if ((j + 1) % self->N == 0) {
            fprintf(f, "\n");
        }       
    }
    fclose(f);

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

void 
sim_compute_wave_stats(sim_t *self, gsl_histogram *density, double *centre, 
        double *width)
{
    unsigned int j;
    unsigned int n = gsl_histogram_bins(density);
    double lower, upper, area, c, w, p, bin;
    gsl_vector *v = gsl_vector_alloc(n);
    for (j = 0; j < n; j++) {
        gsl_histogram_get_range(density, j, &lower, &upper);
        area = M_PI * (upper * upper - lower * lower);
        gsl_vector_set(v, j, gsl_histogram_get(density, j) / area);
    }
    bin = upper - lower;
    /* normalise by the maximum value */
    gsl_vector_scale(v, 1.0 / gsl_vector_max(v));
    c = 0.0;
    w = 0.0;
    for (j = 0; j < n; j++) {
        p = gsl_vector_get(v, j);
        c += p; 
        w += 4 * p * (1 - p);
    }
    gsl_vector_free(v);
    *centre = c * bin;
    *width = w * bin;
}

void
sim_print_short_status(sim_t *self, unsigned long events, gsl_histogram *density,
        gsl_histogram *ancestry) 
{
    /* get the occupancy statistics */
    unsigned int occupied_pixels = 0;
    unsigned int max_occupancy = 0;
    unsigned long total_S_size = 0;
    unsigned long total_C_size = 0;
    unsigned long total_G_size = 0;
    unsigned long total_E_size = 0;
    unsigned int j, o;
    size_t k;
    double w, c, a, sigma;
    
    
    double sim_get_location_variance(sim_t *self, double *z);
    
    for (j = 0; j < gsl_pow_2(self->N); j++) {
        o = avl_count(&self->P[j]);
        if (o > max_occupancy) {
            max_occupancy = o;
        }
        if (o > 0) {
            occupied_pixels++;
        }
        total_S_size += self->total_S_size[j];
        total_C_size += self->total_C_size[j];
        total_G_size += self->total_G_size[j];
        total_E_size += self->total_E_size[j];
        self->total_S_size[j] = 0;
        self->total_C_size[j] = 0;
        self->total_G_size[j] = 0;
        self->total_H_size[j] = 0;
        self->total_E_size[j] = 0;
    }
    w = 0.0;
    c = 0.0;
    sim_compute_wave_stats(self, density, &c, &w);
    gsl_histogram_find(density, c, &k);
    a = gsl_histogram_get(ancestry, k) / gsl_histogram_get(density, k);
    sigma = sim_get_location_variance(self, self->X); 
    assert(total_E_size == self->num_generated_events);
    printf("%ld\t%G\t%d\t%d\t%ld\t%ld\t%ld\t%ld\t%ld\t%d\t%d\t%f\t%f\t%f\t%f\n", 
            events, self->t, self->sample_size, 
            self->ancestral_material, self->num_jumps, 
            self->num_generated_events, total_S_size, 
            total_C_size, total_G_size,
            occupied_pixels, max_occupancy, w, c, a, sigma);
    /* reset the counters */
    self->num_jumps = 0;
    self->num_generated_events = 0;
}


void 
sim_print_individual(sim_t *self, individual_t *ind)
{
    unsigned int last_locus, block_length, ancestral_material;
    avl_node_t *node;
    int_map_value_t *locus_mapping;
    ancestral_material = avl_count(&ind->ancestry);
    printf("%p; am=%d\n", ind, ancestral_material); 
    node = ind->ancestry.head;
    assert(node != NULL);
    locus_mapping = (int_map_value_t *) node->item;
    last_locus = locus_mapping->key;
    block_length = 1;
    while (node->next != NULL) {
        node = node->next; 
        assert(node != NULL);
        locus_mapping = (int_map_value_t *) node->item;
        if (locus_mapping->key == last_locus + 1) {
            block_length++;
        } else {
            printf("\t%d\t%d\n", block_length, last_locus);
            block_length = 1;
        }
        last_locus = locus_mapping->key;
    }
    printf("\t%d\t%d\n", block_length, last_locus);
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

/*
 * Returns the variance in the distance of all individuals 
 * from the point z.
 */
double 
sim_get_location_variance(sim_t *self, double *z)
{
    avl_tree_t chi;
    double d, var;
    unsigned int j;
    uint64_t id, pixel;
    avl_node_t *node, *chi_node;
    individual_t *ind;
    avl_init_tree(&chi, avl_set_compare, NULL); 
    double *v = xmalloc(self->sample_size * sizeof(double)); 
    node = NULL; 
    /* fill chi first */
    for (pixel = 0; pixel < gsl_pow_2(self->N); pixel++) {
        for (node = self->P[pixel].head; node != NULL; node = node->next) {
            id = *((uint64_t *) node->item);
            chi_node = sim_alloc_set_node(self, id); 
            if (avl_insert_node(&chi, chi_node) == NULL) {
                sim_free_avl_set_node(self, chi_node);
            }
        }
    }
    assert(avl_count(&chi) == self->sample_size);
    j = 0;
    for (chi_node = chi.head; chi_node != NULL; chi_node = chi_node->next) {
        id = *((uint64_t *) chi_node->item);
        ind = (individual_t *) id;
        d = sqrt(torus_squared_distance(ind->location, z, self->L));
        v[j] = d;
        j++;
    }
    for (chi_node = chi.head; chi_node != NULL; chi_node = chi_node->next) {
        sim_free_avl_set_node(self, chi_node);
    }
    var = gsl_stats_variance(v, 1, self->sample_size);
    free(v);
    return var;

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
            chi_node = sim_alloc_set_node(self, id); 
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
        for (l = 0; l < self->m; l++) {
            printf("\t%d (%d)\tpi = ", l, self->eta[l]);
            for (j = 0; j < 2 * self->n; j++) {
                printf("%d ", self->pi[l][j]);
            }
            printf("\ttau = ");
            for (j = 0; j < 2 * self->n; j++) {
                printf("%f ", self->tau[l][j]);
            }
            printf("\n");
        }
    }
    for (chi_node = chi.head; chi_node != NULL; chi_node = chi_node->next) {
        sim_free_avl_set_node(self, chi_node);
    }
}

unsigned int 
sim_get_num_blocks(sim_t *self, individual_t *ind)
{
    
    unsigned int last_locus, block_length, ancestral_material, num_blocks;
    avl_node_t *node;
    int_map_value_t *locus_mapping;
    ancestral_material = avl_count(&ind->ancestry);
    printf("%p; am=%d\n", ind, ancestral_material); 
    node = ind->ancestry.head;
    assert(node != NULL);
    locus_mapping = (int_map_value_t *) node->item;
    last_locus = locus_mapping->key;
    block_length = 1;
    num_blocks = 1;
    while (node->next != NULL) {
        node = node->next; 
        assert(node != NULL);
        locus_mapping = (int_map_value_t *) node->item;
        if (locus_mapping->key == last_locus + 1) {
            block_length++;
        } else {
            printf("\t%d\t%d\n", block_length, last_locus);
            block_length = 1;
            num_blocks++;
        }
        last_locus = locus_mapping->key;
    }
    printf("\t%d\t%d\n", block_length, last_locus);
    return num_blocks;
}

void
sim_populate_histograms(sim_t *self, gsl_histogram *density,
        gsl_histogram *ancestry, double *z)
{
    avl_tree_t chi;
    double d;
    unsigned int ancestral_material;
    uint64_t id, pixel;
    avl_node_t *node, *chi_node;
    avl_init_tree(&chi, avl_set_compare, NULL); 
    individual_t *ind;
    node = NULL; 
    /* fill chi first */
    for (pixel = 0; pixel < gsl_pow_2(self->N); pixel++) {
        for (node = self->P[pixel].head; node != NULL; node = node->next) {
            id = *((uint64_t *) node->item);
            chi_node = sim_alloc_set_node(self, id); 
            if (avl_insert_node(&chi, chi_node) == NULL) {
                sim_free_avl_set_node(self, chi_node);
            }
        }
    }
    assert(avl_count(&chi) == self->sample_size);
    gsl_histogram_reset(density);
    gsl_histogram_reset(ancestry);
    for (chi_node = chi.head; chi_node != NULL; chi_node = chi_node->next) {
        id = *((uint64_t *) chi_node->item);
        ind = (individual_t *) id;
        d = sqrt(torus_squared_distance(ind->location, z, self->L));
        if (gsl_histogram_increment(density, d) != 0) {
            //fatal_error("density histogram out of bounnds");
        }
        ancestral_material = avl_count(&ind->ancestry);
        if (gsl_histogram_accumulate(ancestry, d, ancestral_material) != 0) {
            //fatal_error("ancestr histogram out of bounnds");
        }
    }
    for (chi_node = chi.head; chi_node != NULL; chi_node = chi_node->next) {
        sim_free_avl_set_node(self, chi_node);
    }
}


void 
sim_update_hit_distribution(sim_t *self, double *z)
{
    double d = sqrt(torus_squared_distance(self->X, z, self->L));
    if (gsl_histogram_increment(self->hit_rate_distribution, d) != 0) {
        //fatal_error("hit histogram out of bounnds");
    }
}


void 
sim_write_histogram(sim_t *self, gsl_histogram *hist, char *filename)
{
    double lower, upper, area;
    unsigned int j;
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        fatal_error("IO error opening file");
    }
    for (j = 0; j < gsl_histogram_bins(hist); j++) {
        gsl_histogram_get_range(hist, j, &lower, &upper);
        area = M_PI * (upper * upper - lower * lower);
        fprintf(f, "%f\t%f\t%f\n", upper, lower, gsl_histogram_get(hist, j) / area);
    }
    fclose(f);
}

void
sim_write_sample(sim_t *self, char *filename)
{
    avl_tree_t chi;
    double *x;
    uint64_t id, pixel;
    avl_node_t *node, *chi_node;
    avl_init_tree(&chi, avl_set_compare, NULL); 
    individual_t *ind;
    unsigned int first_locus, last_locus, block_length, ancestral_material, 
            num_blocks;
    int_map_value_t *locus_mapping;
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        fatal_error("IO error opening file");
    }
    node = NULL; 
    /* fill chi first */
    for (pixel = 0; pixel < gsl_pow_2(self->N); pixel++) {
        for (node = self->P[pixel].head; node != NULL; node = node->next) {
            id = *((uint64_t *) node->item);
            chi_node = sim_alloc_set_node(self, id); 
            if (avl_insert_node(&chi, chi_node) == NULL) {
                sim_free_avl_set_node(self, chi_node);
            }
        }
    }
    assert(avl_count(&chi) == self->sample_size);
    fprintf(f, "x1\tx2\tancestral_material\tblocks\tfirst\tlast\n");
    for (chi_node = chi.head; chi_node != NULL; chi_node = chi_node->next) {
        id = *((uint64_t *) chi_node->item);
        ind = (individual_t *) id;
        x = ind->location;
        ancestral_material = avl_count(&ind->ancestry);
        node = ind->ancestry.head;
        assert(node != NULL);
        locus_mapping = (int_map_value_t *) node->item;
        last_locus = locus_mapping->key;
        first_locus = last_locus;
        block_length = 1;
        num_blocks = 1;
        while (node->next != NULL) {
            node = node->next; 
            assert(node != NULL);
            locus_mapping = (int_map_value_t *) node->item;
            if (locus_mapping->key == last_locus + 1) {
                block_length++;
            } else {
                block_length = 1;
                num_blocks++;
            }
            last_locus = locus_mapping->key;
        }
        fprintf(f, "%f\t%f\t%d\t%d\t%d\t%d\n", x[0], x[1], 
                ancestral_material, num_blocks, first_locus, last_locus);
    }
    
    for (chi_node = chi.head; chi_node != NULL; chi_node = chi_node->next) {
        sim_free_avl_set_node(self, chi_node);
    }
    fclose(f);
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
#if 0
    assert(self->nu == 2);
    unsigned int k = gsl_ran_binomial(self->rng, self->rho, gap);
    if (k % 2 != 0) {
        ret = (ret + 1) % 2;
    }
#else
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
#endif
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
            self->total_E_size[pixel]++;
            self->num_generated_events++;
        } while (gsl_rng_uniform(self->rng) < 1.0 - jump_proba);
        C_size = gsl_ran_discrete(self->rng, self->beta_distributions[S_size]);
        gsl_ran_choose(self->rng, C, C_size, S, S_size, 
                sizeof(individual_t *));
        //printf("|C| = %d\n", C_size);
        /* update counters */
        sim_update_hit_distribution(self, z);
        self->num_jumps++;
        self->total_H_size[pixel]++;
        self->total_S_size[pixel] += S_size;
        self->total_C_size[pixel] += C_size;
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


/*=============================================================================
 *
 * Algorithm N implementation 
 *
 *=============================================================================
 */



static int 
avl_intset_compare(const void *pa, const void *pb)
{
    int a = *((int *) pa);
    int b = *((int *) pb);
    return (a > b) - (a < b);
}

static void
avl_intset_set_value(avl_node_t *node, int value)
{
    *((int *) node->item) = value;
}

static int  
avl_intset_get_value(avl_node_t *node)
{
    return *((int *) node->item);
}




void 
alg_N_sim_print_state(alg_N_sim_t *self, int extra)
{
    avl_node_t *node;
    unsigned int j, k, c, pixel;
    unsigned int occupied = 0;
    unsigned int coord[2];
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
            coord[0]= j;
            coord[1] = k;
            pixel = pixel_coord_to_index(j, k, self->N); 
            c = avl_count(&self->P[pixel]);
            if (c != 0) {
                occupied++;
            }
            if (c == 0) {
                printf("   ");
            } else {
                printf("%02d ", c);
            }
        }
        printf("|\n");
    }
    printf("    ");
    for (k = 0; k < self->N; k++) {
        printf("___");
    }
    printf("\n");
    if (extra) {
        printf("max_occupancy = %d\n", self->max_occupancy);
        for (j = 0; j <= self->n; j++) {
            printf("occupancy: %d = %d\n", j, avl_count(&self->Q[j]));
            for (k = 0; k < avl_count(&self->Q[j]); k++) {
                node = avl_at(&self->Q[j], k);
                assert(node != NULL);
                pixel = *((int *) node->item);
                index_to_pixel_coord(pixel, self->N, coord);
                printf("\t%d (%d, %d) -> ", pixel, coord[0], coord[1]);
                assert(avl_count(&self->P[pixel]) == j);
                for (node = self->P[pixel].head; node != NULL; 
                        node = node->next) {
                    printf("%d, ", avl_intset_get_value(node)); 
                }
                printf("\n");
            }
        }
        printf("pi  = ");
        for (j = 0; j < 2 * self->n; j++) {
            printf(" %d", self->pi[j]);
        }
        printf("\ntau = ");
        for (j = 0; j < 2 * self->n; j++) {
            printf(" %f", self->tau[j]);
        }
        printf("\n");
    }
}


static void
alg_N_sim_alloc(alg_N_sim_t *self)
{
    unsigned int j, num_nodes;
    self->chi = xmalloc(4 * (self->n + 1) * sizeof(double));
    self->pi = xmalloc(2 * self->n * sizeof(int));
    self->tau = xmalloc(2 * self->n * sizeof(double));
    self->N = self->L / self->s;
    if (fmod(self->L, self->s) != 0.0) {
        fatal_error("L must be a multiple of s");
    }
    self->P = xmalloc(gsl_pow_2(self->N) * sizeof(avl_tree_t));
    for (j = 0; j < gsl_pow_2(self->N); j++) {
        avl_init_tree(&self->P[j], avl_intset_compare, NULL);
    }
    self->rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(self->rng, self->random_seed);
    /* allocate the list memory - each list of length n + 1 */
    /* TODO check if this always works */
    self->max_disc_pixels = 9;
    if (self->s < 1) {
        self->max_disc_pixels = gsl_pow_2(4 * ((int)(self->r / self->s) + 1));
    }
    num_nodes = 2 * self->n * self->max_disc_pixels;
    /* alloc avl node memory */
    self->avl_intset_value_mem = xmalloc(num_nodes * sizeof(int));
    self->avl_intset_node_mem = xmalloc(num_nodes * sizeof(avl_node_t));
    self->avl_intset_node_heap = xmalloc(num_nodes * sizeof(avl_node_t *));
    for (j = 0; j < num_nodes; j++) {
        self->avl_intset_node_heap[j] = self->avl_intset_node_mem + j;
        avl_init_node(self->avl_intset_node_heap[j], 
                &self->avl_intset_value_mem[j]);
    }
    self->avl_intset_node_heap_top = num_nodes - 1; 
    self->Q = xmalloc((self->n + 1) * sizeof(avl_tree_t));
    for (j = 0; j < self->n + 1; j++) {
        avl_init_tree(&self->Q[j], avl_intset_compare, NULL);
    }
}

static void
alg_N_sim_free(alg_N_sim_t *self)
{
    free(self->chi);
    free(self->pi);
    free(self->tau);
    free(self->P);
    free(self->Q);
    free(self->avl_intset_node_heap);
    free(self->avl_intset_node_mem);
    free(self->avl_intset_value_mem);
    gsl_rng_free(self->rng);
}

static void 
alg_N_sim_get_disc_pixels(alg_N_sim_t *self, double *z, unsigned int *p)
{
    get_pixels(self->r, self->s, self->N, z, p);   
}

static avl_node_t * 
alg_N_sim_alloc_avl_intset_node(alg_N_sim_t *self)
{
    avl_node_t *node = NULL;
    if (self->avl_intset_node_heap_top < 0) {
        fatal_error("Out of node space");
    }
    node = self->avl_intset_node_heap[self->avl_intset_node_heap_top];
    self->avl_intset_node_heap_top--;
    return node;

}

static void
alg_N_sim_free_avl_intset_node(alg_N_sim_t *self, avl_node_t *node) 
{
    self->avl_intset_node_heap_top++;
    self->avl_intset_node_heap[self->avl_intset_node_heap_top] = node;
}


static void
alg_N_sim_insert_lineage(alg_N_sim_t *self, unsigned int lineage, 
        unsigned int pixel)
{
    unsigned int h; 
    avl_node_t *node;
    h = avl_count(&self->P[pixel]);
    node = alg_N_sim_alloc_avl_intset_node(self);
    avl_intset_set_value(node, lineage);
    avl_insert_node(&self->P[pixel], node);
    if (h == 0) {
        node = alg_N_sim_alloc_avl_intset_node(self);
    } else {
        node = avl_search(&self->Q[h], &pixel);
        assert(node != NULL);
        avl_unlink_node(&self->Q[h], node); 
    }
    h++;
    avl_intset_set_value(node, pixel);
    avl_insert_node(&self->Q[h], node);
    if (h > self->max_occupancy) {
        //printf("increased max_occupancy: %d\n", h);
        self->max_occupancy =  h;
    }
}

static void
alg_N_sim_remove_lineage(alg_N_sim_t *self, unsigned int lineage, 
        unsigned int pixel)
{
    avl_node_t *node;
    unsigned int h = avl_count(&self->P[pixel]);
    node = avl_search(&self->P[pixel], &lineage);
    assert(node != NULL);
    avl_unlink_node(&self->P[pixel], node); 
    alg_N_sim_free_avl_intset_node(self, node);
    /* update occupancy */
    node = avl_search(&self->Q[h], &pixel);
    assert(node != NULL);
    avl_unlink_node(&self->Q[h], node); 
    if (h == 1) {
        alg_N_sim_free_avl_intset_node(self, node);
    } else {
        avl_intset_set_value(node, pixel);
        avl_insert_node(&self->Q[h - 1], node);
    }
}

void 
alg_N_sim_run(alg_N_sim_t *self)
{
    unsigned int h, j, k, v, parent, kappa, eta;
    unsigned int *pixels = xmalloc((self->max_disc_pixels + 1) * sizeof(int));
    unsigned int *S = xmalloc(self->n * self->max_disc_pixels * sizeof(int));
    unsigned int *C = xmalloc((self->n + 1) * sizeof(int));
    double r2 = gsl_pow_2(self->r);
    double *p = xmalloc(self->n * sizeof(double)); 
    double *b = xmalloc(self->n * sizeof(double)); 
    double *ubar = xmalloc((self->n + 1) * sizeof(double));
    double t, w, Lambda_const, z[2], *x, jump_proba;
    double u = self->u;
    avl_node_t *node;
    unsigned int coord[2];
    for (j = 0; j < 2 * self->n; j++) {
        self->pi[j] = 0;
        self->tau[j] = 0.0;
        self->chi[2 * j] = -1.0;
        self->chi[2 * j + 1] = -1.0;
    }
    for (j = 1; j <= self->n; j++) {
        self->chi[2 * j] = self->X[2 * (j - 1)]; 
        self->chi[2 * j + 1] = self->X[2 * (j - 1) + 1]; 
    }
    for (j = 0; j <= self->n; j++) {
        ubar[j] = 1.0 - gsl_pow_int(1.0 - u, j);
    }
    kappa = self->n;
    eta = self->n + 1;
    t = 0.0;
    self->max_occupancy = 0;
    for (j = 1; j <= self->n; j++) {
        alg_N_sim_get_disc_pixels(self, &self->chi[2 * j], pixels); 
        for (k = 1; k <= pixels[0]; k++) {
            alg_N_sim_insert_lineage(self, j, pixels[k]); 
        }
    }
    Lambda_const = gsl_pow_2(self->s / self->L);
    self->num_jumps = 0;
    self->num_generated_events = 0;
    self->total_max_occupancy = 0;
    while (kappa > 1 && self->num_jumps < self->max_jumps) {
        //alg_N_sim_print_state(self, 2);
        /* adjust max_occupancy downwards, if necessary */
        while (avl_count(&self->Q[self->max_occupancy]) == 0) {
            self->max_occupancy--;
            assert(self->max_occupancy > 0);
            //printf("decreased max_occupancy: %d\n", self->max_occupancy);
        }
        self->total_max_occupancy += self->max_occupancy;
        w = 0.0;
        p[0] = 0.0;
        for (h = 1; h <= self->max_occupancy; h++) {
            p[h] = avl_count(&self->Q[h]) * ubar[h]; 
            w += p[h];
        }
        /* normalise p */
        for (h = 1; h <= self->max_occupancy; h++) {
            p[h] /= w;
        }
        do {
            t += gsl_ran_exponential(self->rng, 1.0 / (w * Lambda_const));
            h = probability_list_select(p, self->max_occupancy + 1, 
                    gsl_rng_uniform(self->rng));
            k = gsl_rng_uniform_int(self->rng, avl_count(&self->Q[h]));
            node = avl_at(&self->Q[h], k); 
            v = avl_intset_get_value(node);
            index_to_pixel_coord(v, self->N, coord);
            z[0] = self->s * (coord[0] + gsl_rng_uniform(self->rng));
            z[1] = self->s * (coord[1] + gsl_rng_uniform(self->rng));
            assert(avl_count(&self->P[v]) != 0);
            S[0] = 0;
            for (node = self->P[v].head; node != NULL; node = node->next) {
                j = avl_intset_get_value(node); 
                x = &self->chi[2 * j];
                if (torus_squared_distance(z, x, self->L) < r2) {
                    S[0]++;
                    S[S[0]] = j; 
                }
            }
            jump_proba = ubar[S[0]] / ubar[h]; 
            assert(jump_proba <= 1.0);
            self->num_generated_events++;
        } while (gsl_rng_uniform(self->rng) < 1.0 - jump_proba);
        self->num_jumps++;
        for (j = 1; j <= S[0]; j++) {
            b[j - 1] = beta(S[0], j, u);
        }
        C[0] = 1 + probability_list_select(b, S[0], gsl_rng_uniform(self->rng));
        gsl_ran_choose(self->rng, C + 1, C[0], S + 1, S[0], sizeof(int));
        parent = C[0] == 1 ? C[1] : eta;
        for (j = 1; j <= C[0]; j++) {
            x = &self->chi[2 * C[j]];
            alg_N_sim_get_disc_pixels(self, x, pixels);
            for (k = 1; k <= pixels[0]; k++) {
                alg_N_sim_remove_lineage(self, C[j], pixels[k]); 
            }
        }
        x = &self->chi[2 * parent];
        random_point_torus_disc(x, z, self->r, self->L, self->rng);
        alg_N_sim_get_disc_pixels(self, x, pixels);
        for (k = 1; k <= pixels[0]; k++) {
            alg_N_sim_insert_lineage(self, parent, pixels[k]);
        }
        if (C[0] > 1) {
            for (j = 1; j <= C[0]; j++) {
                self->pi[C[j]] = eta;
            }
            self->tau[eta] = t;
            eta = eta + 1;
            kappa = kappa - C[0] + 1;
        }
    }
    free(pixels);
    free(S);
    free(C);
    free(p);
    free(b);
    free(ubar);
}





/*=============================================================================
 *
 * Run functions. 
 *
 *=============================================================================
 */

void 
run_algorithm_M_identity(double u, double x)
{
    double mu = 1e-10;
    double accuracy_goal = 1e-4;
    double X[] = {0.0, 0.0, 0.0, 0.0};
    unsigned int l;
    unsigned int coalesced_loci = 0;
    unsigned int num_distinct_times = 0;
    double last_t;
    double F = 0.0;
    double *t;
    char filename[1024];
    FILE *f;
    sim_t sim;
    sim.m = (int) 1e6;
    sim.nu = 2;
    sim.n = sizeof(X) / (2 * sizeof(double));
    sim.random_seed = 1;
    sim.r = 1.0;
    sim.u = u;
    sim.lambda = 1.0;
    sim.rho = 1e-2;
    sim.s = 2.5;
    sim.L = 1000; 
    sim.X = X;
    sim.beta_threshold = 1e-14;
    sim.max_time = log(accuracy_goal) / (-2 * mu);
    sim.max_sample_size = (unsigned int) 1e6;
    sim.max_jumps = ULLONG_MAX;
    X[0] = x;
    sprintf(filename, "M-identity/%f_%f.dat", u, x); 
    f = fopen(filename, "w");
    if (f == NULL) {
        fatal_error("cannot open file");
    }
    sim_set_max_occupancy(&sim, 0.5);
    sim_alloc(&sim);
    sim_initialise(&sim);
    sim_run(&sim);
    
    t = xmalloc(sim.m * sizeof(double)); 
    for (l = 0; l < sim.m; l++) {
        if (sim.pi[l][1] != 0) {
            t[coalesced_loci] = sim.tau[l][3];
            coalesced_loci++;
        }
    }
    if (coalesced_loci > 0) {
        gsl_sort(t, 1, coalesced_loci);
        last_t = -1.0; 
        for (l = 0; l < coalesced_loci; l++) {
            if (t[l] != last_t) {
                last_t = t[l];
                F += exp(-2 * mu * t[l]);
                num_distinct_times++;
            }
        }
        F /= num_distinct_times;
    }
    fprintf(f, "%f\t%f\t%.14f\t%d\n", u, x, F, num_distinct_times);
    fclose(f);
    sim_free(&sim);
    free(t);
}

    
void 
run_algorithm_M_stats(double u, long seed)
{
    char filename[1024];
    int j;
    double X[] = {50.5, 50.5, 50.5, 50.5};
    unsigned long events = 0;
    gsl_histogram *density = gsl_histogram_alloc(100);
    gsl_histogram *ancestry = gsl_histogram_alloc(100);
    sim_t sim;
    sim.m = (int) 1e6;
    sim.L = 1000; 
    sim.nu = 2;
    sim.n = sizeof(X) / (2 * sizeof(double));
    sim.random_seed = seed;
    sim.r = 1.0;
    sim.u = u;
    sim.lambda = 1.0;
    sim.rho = 0.5; 
    sim.s = 1.0;
    sim.X = X;
    sim.beta_threshold = 1e-14;
    sim.max_time = DBL_MAX;
    sim.max_sample_size = (unsigned int) 2e6;
    sim_set_max_occupancy(&sim, 0.5);
    sim_alloc(&sim);
    sim_initialise(&sim);
    sim.max_jumps = (unsigned long) 1e2;
    sim_print_parameters(&sim);
    gsl_histogram_set_ranges_uniform(density, 0.0, 100.0);
    gsl_histogram_set_ranges_uniform(ancestry, 0.0, 100.0);
    printf("t   T   n   kappa   jumps   events  S   C   G    "
            "occupied_pixels max_occupancy w c Ac sigma\n");
    for (j = 1; j <  100; j++) {
        sim_run(&sim);
        events += sim.max_jumps;
        //sim_print_state(&sim, 0);
        /*
        sprintf(filename, "output_NOBACKUP_/H_%ld.dat", events);
        sim_write_distribution(&sim, filename, sim.total_H_size);
        sprintf(filename, "output_NOBACKUP_/E_%ld.dat", events);
        sim_write_distribution(&sim, filename, sim.total_E_size);
        sprintf(filename, "output_NOBACKUP_/C_%ld.dat", events);
        sim_write_distribution(&sim, filename, sim.total_C_size);
        sprintf(filename, "output_NOBACKUP_/S_%ld.dat", events);
        sim_write_distribution(&sim, filename, sim.total_S_size);
        sprintf(filename, "output_NOBACKUP_/G_%ld.dat", events);
        sim_write_distribution(&sim, filename, sim.total_G_size);
        */
        sim_populate_histograms(&sim, density, ancestry, X);
        sim_print_short_status(&sim, events, density, ancestry);
        /*
        sprintf(filename, "output_NOBACKUP_/occupancy_%ld.dat", events);
        sim_write_occupancy_distribution(&sim, filename);
        sprintf(filename, "output_NOBACKUP_/ancestry_%ld.dat", events);
        sim_write_ancestry_distribution(&sim, filename);
        sprintf(filename, "output_NOBACKUP_/blocks_%ld.dat", events);
        sim_write_block_distribution(&sim, filename);
        sprintf(filename, "output_NOBACKUP_/sample_%ld.dat", events);
        sim_write_sample(&sim, filename);
        */
        sprintf(filename, "output_NOBACKUP_/density_hist_%ld.dat", events);
        sim_write_histogram(&sim, density, filename);
        sprintf(filename, "output_NOBACKUP_/ancestry_hist_%ld.dat", events);
        sim_write_histogram(&sim, ancestry, filename);
        /*
        sprintf(filename, "output_NOBACKUP_/hitrate_hist_%ld.dat", events);
        sim_write_histogram(&sim, sim.hit_rate_distribution, filename);
        */
        gsl_histogram_reset(sim.hit_rate_distribution);
        //sprintf(filename, "output_NOBACKUP_/occupancy_%d.png", j);
        //sim_draw_occupancy(&sim, filename);
        //sprintf(filename, "output_NOBACKUP_/ancestral_material_%d.png", j);
        //sim_draw_ancestral_material(&sim, filename);
    }
    sim_free(&sim);
}

void 
run_algorithm_N_identity(double u, long seed)
{
    int mrca;
    double t, F;
    unsigned int j;
    unsigned int n = 100;
    double max_x = 10; 
    double mu = 1e-6;
    double accuracy_goal = 1e-6;
    double *X = xmalloc(2 * n * sizeof(double));
    char filename[8192];
    FILE *f;
    sv_t sv;
    alg_N_sim_t sim;
    sim.n = n; 
    sim.random_seed = seed;
    sim.r = 1.0;
    sim.u = u;
    sim.lambda = 1.0;
    sim.s = 2;
    sim.L = 100; 
    sim.X = X;
    sim.max_time = log(accuracy_goal) / (-2 * sim.mu);
    sim.max_jumps = ULLONG_MAX;
    sprintf(filename, "output_NOBACKUP_/%.8f_%ld.dat", u, seed);
    f = fopen(filename, "w"); 
    if (f == NULL) {
        fatal_error("Cannot open file");
    }
    
    alg_N_sim_alloc(&sim);
    
    X[0] = 0.0;
    X[1] = 0.0;
    X[2] = 0.0;
    X[3] = 0.0;
    for (j = 2; j < n; j++) {
        X[2 * j] = 0.0;
        X[2 * j + 1] = (j - 1) * max_x / n;
    }
    alg_N_sim_run(&sim);
    sv_init(&sv, 2 * n, sim.pi); 
    for (j = 1; j < n; j++) {
        mrca = sv_nca(&sv, 1, j + 1);
        t = sim.tau[mrca];
        F = exp(-2 * mu * t);
        fprintf(f, "%ld\t%f\t%.8f\n", seed, X[2 * j + 1], F);
    }
    sv_free(&sv);
    alg_N_sim_free(&sim);
    free(X);
    fclose(f);
}


void 
run_algorithm_N_analysis(int n, unsigned long max_jumps, long seed)
{
    int j;
    double *X = malloc(2 * n * sizeof(double)); 
    double s, user_time;
    struct tms before, after;
    alg_N_sim_t sim;
    sim.n = n; 
    sim.random_seed = seed;
    sim.r = 1.0;
    sim.u = 0.00125;
    sim.lambda = 1.0;
    sim.mu = 1e-10;
    sim.s = 2.0;
    sim.X = X;
    sim.max_time = DBL_MAX; 
    sim.max_jumps = max_jumps;
    printf("n\ts\tj\tt\tg\to\n"); 
    for (s = 1.0; s < 4.25; s += 0.125) { 
        sim.s = s;
        sim.L = 10000 * s; 
        alg_N_sim_alloc(&sim);
        for (j = 0; j < 2 * n; j++) {
            X[j] = sim.L * gsl_rng_uniform(sim.rng);
        }
        times(&before);
        alg_N_sim_run(&sim);
        times(&after);
        user_time = ((double) (after.tms_utime - before.tms_utime)) / sysconf(_SC_CLK_TCK);
        printf("%d\t%f\t%ld\t%.8f", n, s, sim.num_jumps, user_time);
        printf("\t%ld\t%ld\n", sim.num_generated_events, sim.total_max_occupancy);
        alg_N_sim_free(&sim);
    }
    free(X);
}

void 
run_algorithm_M_analysis(double u, unsigned long max_jumps, long seed)
{
    unsigned int j;
    double *X = malloc(4 * sizeof(double)); 
    double s, user_time, N;
    struct tms before, after;
    sim_t sim;
    sim.m = (unsigned int) 1e6; 
    sim.nu = 2;
    sim.n = 2; 
    sim.random_seed = seed;
    sim.r = 1.0;
    sim.u = u;
    sim.lambda = 1.0;
    sim.rho = 0.5; 
    sim.s = 0.1;
    sim.X = X;
    sim.beta_threshold = 1e-8;
    sim.max_time = DBL_MAX; 
    sim.max_sample_size = (unsigned int) 1e6;
    sim.max_jumps = max_jumps;
    N = sim.nu / sim.u;
    printf("# N = %f\n", N);
    printf("s\tjumps\ttime\tsample_size\n");
    for (s = 0.5; s <= 2.0; s += 0.125) { 
        sim.s = s;
        sim_set_max_occupancy(&sim, 0.5);
        sim.L = 1000 * s; 
        sim_alloc(&sim);
        for (j = 0; j < 2 * sim.n; j++) {
            X[j] = sim.L * gsl_rng_uniform(sim.rng);
        }
        sim_initialise(&sim);
        times(&before);
        sim_run(&sim);
        times(&after);
        user_time = ((double) (after.tms_utime - before.tms_utime)) / sysconf(_SC_CLK_TCK);
        printf("%f\t%ld\t%.8f\t%d\n", s, sim.num_jumps, user_time, sim.sample_size);
        sim_free(&sim);
    }
}




int 
main(int argc, char **argv)
{
    char *cmd;
    double u, x;
    long seed;
    int n;
    unsigned long max_jumps;
    if (argc < 2) {
        fatal_error("arguments required");
    }
    cmd = argv[1];
    /*
    if (strstr(cmd, "identity") != NULL) { 
       //run_identity(num_replicates);
    } else if (strstr(cmd, "analysis") != NULL) { 
        //run_analysis();
    } else if (strstr(cmd, "pairwise") != NULL) { 
        //run_pairwise(num_replicates);
    */
    if (strstr(cmd, "N-identity") != NULL) { 
        if (argc != 4) {
            fatal_error("usage: N-identity u seed");
        }
        u = atof(argv[2]);
        seed = atoi(argv[3]);
        run_algorithm_N_identity(u, seed);
    } else if (strstr(cmd, "N-analysis") != NULL) {
        if (argc != 5) {
            fatal_error("usage: N-analysis sample_size max_jumps seed");
        }
        n = atoi(argv[2]);
        max_jumps = (unsigned long) atof(argv[3]);
        seed = atoi(argv[4]);
        run_algorithm_N_analysis(n, max_jumps, seed);
    } else if (strstr(cmd, "M-stats") != NULL) { 
        if (argc != 4) {
            fatal_error("usage: M-stats u seed");
        }
        u = atof(argv[2]);
        seed = atoi(argv[3]);
        run_algorithm_M_stats(u, seed);
    } else if (strstr(cmd, "M-identity") != NULL) { 
        if (argc != 4) {
            fatal_error("usage: M-identity u x");
        }
        u = atof(argv[2]);
        x = atof(argv[3]);
        run_algorithm_M_identity(u, x);
     } else if (strstr(cmd, "M-analysis") != NULL) {
        if (argc != 5) {
            fatal_error("usage: M-analysis u max_jumps seed");
        }
        u = atof(argv[2]);
        max_jumps = (unsigned long) atof(argv[3]);
        seed = atoi(argv[4]);
        run_algorithm_M_analysis(u, max_jumps, seed);
    } else {
        fatal_error("unrecognised argument");
    }
    return 0;
}

