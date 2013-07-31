/*
** Copyright (C) 2013 Jerome Kelleher <jerome.kelleher@ed.ac.uk>
**  
** This file is part of sim.
** 
** sim is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
** 
** sim is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with sim.  If not, see <http://www.gnu.org/licenses/>.
*/

/* 
 * Simple command line front end for the sim code. This is intended 
 * for development purposes only!
 */

#include "sim.h"
#include "util.h"

#include <limits.h>
#include <string.h>
#include <stdarg.h>
#include <libconfig.h>


/* From version 1.4, libconfig has swapped long for int as the integer type.
 * The LIBCONFIG_VER_MAJOR macro was introduced at version 1.4 also, so this 
 * should be a safe test.
 */
#ifdef LIBCONFIG_VER_MAJOR
typedef int libconfig_int;
#else
typedef long libconfig_int;
#endif


static void 
fatal_error(const char *msg, ...)
{
    va_list argp;
    fprintf(stderr, "sim:");
    va_start(argp, msg);
    vfprintf(stderr, msg, argp);
    va_end(argp);
    fprintf(stderr, "\n");
    exit(EXIT_FAILURE);
}
/*
static void
read_recombination_probabilities(sim_t *self, config_t *config)
{
    int j;
    int m;
    const char *key = "recombination_probabilities";
    config_setting_t *s;
    config_setting_t *setting = config_lookup(config, key); 
    
    if (setting == NULL) {
        fatal_error(key, " is a required parameter");
    }
    if (config_setting_is_array(setting) == CONFIG_FALSE) {
        fatal_error(key, "must be an array");
    }
    m = config_setting_length(setting) + 1;
    self->num_loci = m;
    self->recombination_probabilities = xmalloc(m * sizeof(double));
    for (j = 0; j < m - 1; j++) {
        s = config_setting_get_elem(setting, j);
        if (s == NULL) {
            fatal_error("error reading %s[%d]", key, j);
        }
        self->recombination_probabilities[j] = config_setting_get_float(s);
    }
    self->recombination_probabilities[m - 1] = 0.0;
}
*/


/*
static void
read_sample(sim_t *self, config_t *config)
{
    int j;
    config_setting_t *s;
    config_setting_t *setting = config_lookup(config, "sample"); 
    if (setting == NULL) {
        fatal_error("sample is a required parameter");
    }
    if (config_setting_is_list(setting) == CONFIG_FALSE) {
        fatal_error("sample must be a list");
    }
    self->sample_size = config_setting_length(setting);
    if (self->sample_size < 1) {
        fatal_error("sample size must be > 0");
    }
    self->sample = xmalloc(2 * self->sample_size * sizeof(double));
    for (j = 0; j < self->sample_size; j++) {
        s = config_setting_get_elem(setting, j);
        if (s == NULL) {
            fatal_error("error reading sample[%d]", j);
        }
        if (config_setting_is_array(s) == CONFIG_FALSE) {
            fatal_error("sample[%d] not an array", j);
        }
        if (config_setting_length(s) != 2) { 
            fatal_error("sample[%d] not 2D", j);
        }
        self->sample[2 * j] = config_setting_get_float_elem(s, 0);
        self->sample[2 * j + 1] = config_setting_get_float_elem(s, 1);
    }
}
*/
/*
static void
read_events(sim_t *self, config_t *config)
{
    double u, r, rate, theta, alpha; 
    int j;
    int e;
    const char *type;
    config_setting_t *s, *t;
    config_setting_t *setting = config_lookup(config, "events"); 
    if (setting == NULL) {
        fatal_error("events is a required parameter");
    }
    if (config_setting_is_list(setting) == CONFIG_FALSE) {
        fatal_error("events must be a list");
    }
    self->num_event_classes = config_setting_length(setting);
    e = self->num_event_classes;
    if (e < 1) {
        fatal_error("events must be > 0");
    }
    self->event_classes = xmalloc(e * sizeof(event_class_t));
    for (j = 0; j < (int) e; j++) {
        s = config_setting_get_elem(setting, j);
        if (s == NULL) {
            fatal_error("error reading events[%d]", j);
        }
        if (config_setting_is_group(s) == CONFIG_FALSE) {
            fatal_error("events[%d] not a group", j);
        }
        t = config_setting_get_member(s, "rate"); 
        if (t == NULL) {
            fatal_error("rate not specified");
        }
        rate = config_setting_get_float(t);
        if (rate < 0.0) {
            fatal_error("event rate must be > 0");
        }
        t = config_setting_get_member(s, "type"); 
        if (t == NULL) {
            fatal_error("type not specified");
        }
        type = config_setting_get_string(t);
        if (strcmp(type, "disc") == 0) {
            t = config_setting_get_member(s, "u"); 
            if (t == NULL) {
                fatal_error("u not specified");
            }
            u = config_setting_get_float(t);
            t = config_setting_get_member(s, "r"); 
            if (t == NULL) {
                fatal_error("r not specified");
            }
            r = config_setting_get_float(t);
            alloc_disc_event_class(&self->event_classes[j], rate, r, u);
        } else if (strcmp(type, "gauss") == 0) {
            t = config_setting_get_member(s, "u0"); 
            if (t == NULL) {
                fatal_error("u0 not specified");
            }
            u = config_setting_get_float(t);
            t = config_setting_get_member(s, "theta"); 
            if (t == NULL) {
                fatal_error("theta not specified");
            }
            theta = config_setting_get_float(t);
            t = config_setting_get_member(s, "alpha"); 
            if (t == NULL) {
                fatal_error("alpha not specified");
            }
            alpha = config_setting_get_float(t);
            alloc_gaussian_event_class(&self->event_classes[j], rate, theta,
                    alpha, u);

        } else {
            fatal_error("unknown event type '%s'", type);
        }
        
    }
}

static void 
read_config(sim_t *self, const char *filename)
{
    int err;
    libconfig_int tmp;
    config_t *config = xmalloc(sizeof(config_t)); 
    config_init(config);
    err = config_read_file(config, filename);
    if (err == CONFIG_FALSE) {
        fatal_error("configuration error:%s at line %d in file %s\n", 
                config_error_text(config), config_error_line(config), 
                filename);
    }
    if (config_lookup_int(config, "num_parents", &tmp) 
            == CONFIG_FALSE) {
        fatal_error("num_parents is a required parameter");
    }
    self->num_parents = tmp;
    if (config_lookup_int(config, "random_seed", &tmp) 
            == CONFIG_FALSE) {
        fatal_error("random_seed is a required parameter");
    }
    self->random_seed = tmp;
    if (config_lookup_int(config, "kdtree_bucket_size", &tmp) 
            == CONFIG_FALSE) {
        fatal_error("kdtree_bucket_size is a required parameter");
    }
    self->kdtree_bucket_size = tmp;
    if (config_lookup_int(config, "max_kdtree_insertions", &tmp) 
                == CONFIG_FALSE) {
        fatal_error("max_kdtree_insertions is a required parameter");
    }
    self->max_kdtree_insertions = tmp; 
    if (config_lookup_int(config, "max_lineages", &tmp) 
                == CONFIG_FALSE) {
        fatal_error("max_lineages is a required parameter");
    }
    self->max_lineages = tmp; 
    if (config_lookup_float(config, "torus_diameter", &self->torus_diameter) 
            == CONFIG_FALSE) {
        fatal_error("torus_diameter is a required parameter");
    }
    if (config_lookup_float(config, "max_time", 
                &self->max_time) == CONFIG_FALSE) {
        fatal_error("max_time is a required parameter");
    }
    read_events(self, config);    
    //read_sample(self, config);
    //read_recombination_probabilities(self, config); 

    config_destroy(config);
    free(config);
}

*/


int
main(int argc, char** argv)
{
    //int ret;
    //int not_done = 1;
    sim_t *self = xcalloc(1, sizeof(sim_t));
    printf("%s", argv[0]);  // REMOVE
    if (argc != 2) {
        fatal_error("usage: sim <configuration file>");
    }
    //read_config(self, argv[1]); 
    /*
    ret = sim_initialise(self);
    ERCS_ERROR_CHECK(ret, out); 
    while (not_done) {
        ret = sim_simulate(self, UINT_MAX);
        ERCS_ERROR_CHECK(ret, out);
        not_done = ret == ERCS_SIM_NOT_DONE;
    }
    sim_print_state(self);
   
out:
    if (ret < 0) {
        printf("Error occured: %d: %s\n", ret, sim_error_str(ret));
    }
    */
    sim_free(self);
    /*
    free(self->recombination_probabilities);
    free(self->event_classes);
    free(self->sample); 
    */
    free(self);

    return EXIT_SUCCESS;
}
