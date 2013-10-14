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
#include <Python.h>
#include <structmember.h>

#include <float.h>
#include <gsl/gsl_math.h>
#include "lib/util.h"
#include "lib/sim.h"
#include "lib/nystrom.h"

#if PY_MAJOR_VERSION >= 3
#define IS_PY3K
#endif

#define MODULE_DOC \
"Low level interface for discsim" 

static PyObject *DiscsimInputError;
static PyObject *DiscsimLibraryError;

typedef struct {
    PyObject_HEAD
    sim_t *sim;
} Simulator;

typedef struct {
    PyObject_HEAD
    nystrom_t *nystrom;
} IdentitySolver;

static void 
handle_library_error(int err)
{
    PyErr_SetString(DiscsimLibraryError, discsim_error_message(err));
}

static void 
handle_input_error(const char *err)
{
    PyErr_SetString(DiscsimInputError, err);
}

/*
 * Retrieves a number value with the specified key from the specified 
 * dictionary.
 */
static PyObject *
get_dict_number(PyObject *dict, const char *key_str)
{
    PyObject *ret = NULL;
    PyObject *value;
    PyObject *key = Py_BuildValue("s", key_str);
    if (!PyDict_Contains(dict, key)) {
        PyErr_Format(DiscsimInputError, "'%s' not specified", key_str); 
        goto out;
    }
    value = PyDict_GetItem(dict, key);
    if (!PyNumber_Check(value)) {
        PyErr_Format(DiscsimInputError, "'%s' is not number", key_str); 
        goto out;
    }
    ret = value;
out:
    Py_DECREF(key);
    return ret;
}



static int 
discsim_parse_event_classes(PyObject *py_events, event_class_t *events) 
{
    int ret = -1;
    int j, size;
    double rate, u, r;
    PyObject *item, *value;
    size = PyList_Size(py_events);
    if (size == 0) {
        PyErr_SetString(DiscsimInputError, "must have > 0 events"); 
        goto out;
    }
    for (j = 0; j < size; j++) {
        item = PyList_GetItem(py_events, j);
        if (!PyDict_Check(item)) {
            PyErr_SetString(DiscsimInputError, "not a dictionary"); 
            goto out;
        }
        value = get_dict_number(item, "rate");
        if (value == NULL) {
            goto out;
        }
        rate = PyFloat_AsDouble(value);
        value = get_dict_number(item, "r");
        if (value == NULL) {
            goto out;
        }
        r = PyFloat_AsDouble(value);
        value = get_dict_number(item, "u");
        if (value == NULL) {
            goto out;
        }
        u = PyFloat_AsDouble(value);
        events[j].rate = rate;
        events[j].r = r;
        events[j].u = u;
    }
    ret = 0;
out:
    return ret;
}

/*===================================================================
 * Simulator
 *===================================================================
 */
static int 
Simulator_check_sim(Simulator *self) 
{
    int ret = 0;
    if (self->sim == NULL) {
        PyErr_SetString(PyExc_SystemError, "simulator not initialised");
        ret = -1; 
    }
    return ret;
}


static int
Simulator_parse_sample(Simulator *self, PyObject *py_sample)
{
    int ret = -1;
    int size;
    int j, k;
    double v;
    PyObject *item, *value;    
    int n = PyList_Size(py_sample);
    if (n == 0) {
        PyErr_SetString(DiscsimInputError, "Empty sample"); 
        goto out;
    }
    self->sim->sample_size = n; 
    self->sim->sample = PyMem_Malloc(2 * n * sizeof(double));
    memset(self->sim->sample, 0, 2 * n * sizeof(double));
    for (j = 0; j < n; j++) {
        item = PyList_GetItem(py_sample, j);
        if (self->sim->dimension == 1) {
            value = item;
            if (!PyNumber_Check(value)) {
                PyErr_SetString(DiscsimInputError, 
                        "Locations must be numeric");
                goto out;
            }
            v = PyFloat_AsDouble(value);
            self->sim->sample[j * 2] = v;
            self->sim->sample[j * 2 + 1] = 0.0;
            if (v < 0.0 || v >= self->sim->torus_diameter) {
                PyErr_SetString(DiscsimInputError, 
                        "sample location: must have 0 <= v < L"); 
                goto out;
            }
        } else {
            size = 0;
            if (!PyTuple_Check(item)) {
                PyErr_SetString(DiscsimInputError, "Samples must be 2-tuples"); 
                goto out;
            } else {
                size = PyTuple_Size(item);
                if (size != 2) {
                    PyErr_SetString(DiscsimInputError, "Dimension != 2 not supported"); 
                    goto out;
                }
                for (k = 0; k < 2; k++) {
                    value = PyTuple_GetItem(item, k);
                    if (!PyNumber_Check(value)) {
                        PyErr_SetString(DiscsimInputError, 
                                "Locations must be numeric");
                        goto out;
                    }
                    v = PyFloat_AsDouble(value);
                    self->sim->sample[j * 2 + k] = v;
                    if (v < 0.0 || v >= self->sim->torus_diameter) {
                        PyErr_SetString(DiscsimInputError, 
                                "sample location: must have 0 <= v < L"); 
                        goto out;
                    }
                }
            }
        }
    }
    ret = 0;
out:   
    return ret; 
}

static int
Simulator_parse_events(Simulator *self, PyObject *py_events)
{
    int ret = -1;
    int size;
    size = PyList_Size(py_events);
    if (size == 0) {
        PyErr_SetString(DiscsimInputError, "must have > 0 events"); 
        goto out;
    }
    self->sim->num_event_classes = size; 
    self->sim->event_classes = PyMem_Malloc(size * sizeof(event_class_t));
    if (self->sim->event_classes == NULL) {
        ret = ERR_ALLOC_FAILED; 
        goto out;
    }
    ret = discsim_parse_event_classes(py_events, self->sim->event_classes);
out:
    return ret; 
}
   
    
static void
Simulator_dealloc(Simulator* self)
{
    if (self->sim != NULL) {
        if (self->sim->sample != NULL) {
            PyMem_Free(self->sim->sample);
        }
        if (self->sim->event_classes != NULL) {
            PyMem_Free(self->sim->event_classes);
        }
        sim_free(self->sim);
        PyMem_Free(self->sim);
        self->sim = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int 
Simulator_check_input(Simulator *self)
{
    int ret = -1;
    unsigned int j;
    sim_t *sim = self->sim;
    event_class_t *e;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    if (sim->dimension < 1 || sim->dimension > 2) {
        handle_input_error("dimension must be 1 or 2");
        goto out;
    }
    if (sim->dimension == 1 && sim->pixel_size != 2.0) {
        handle_input_error("pixel size must be 2.0 in 1D");
        goto out;
    }
    if (sim->simulate_pedigree < 0 || sim->simulate_pedigree > 1) {
        handle_input_error("simulate_pedigree must be 0 or 1");
        goto out;
    }
    if (sim->simulate_pedigree == 1 && sim->num_loci != 1) {
        handle_input_error("m must be 1 for pedigree simulation");
        goto out;
    }
    if (sim->torus_diameter <= 0.0) {
        handle_input_error("must have torus_edge > 0");
        goto out;
    }
    if (sim->num_loci == 0) {
        handle_input_error("must have num_loci > 0");
        goto out;
    }
    if (sim->num_parents == 0) {
        handle_input_error("must have num_parents > 0");
        goto out;
    }
    if (sim->max_population_size == 0) {
        handle_input_error("must have max_population_size > 0");
        goto out;
    }
    if (sim->max_occupancy == 0) {
        handle_input_error("must have max_occupancy > 0");
        goto out;
    }
    if (sim->recombination_probability < 0 || 
            sim->recombination_probability > 1) {
        handle_input_error("must have 0 <= recombination_probability <= 1");
        goto out;
    }
    if (sim->pixel_size <= 0 || sim->pixel_size > sim->torus_diameter / 4) {
        handle_input_error("must have 0 < pixel_size <= L/4 ");
        goto out;
    }
    if (fmod(sim->torus_diameter, sim->pixel_size) != 0.0) {
        handle_input_error("L/s must be an integer");
        goto out;
    }
    if (sim->num_event_classes == 0) {
        handle_input_error("at least one event class required");
        goto out;
    }
    for (j = 0; j < sim->num_event_classes; j++) {
        e = &sim->event_classes[j];
        if (e->r <= 0.0 || e->r > sim->torus_diameter / 4.0) {
            handle_input_error("must have 0 < r < L / 4");
            goto out;
        }
        if (e->u <= 0.0 || e->u >= 1.0) {
            handle_input_error("must have 0 < u < 1");
            goto out;
        }
        if (e->rate <= 0.0) {
            handle_input_error("must have 0 < rate < 1");
            goto out;
        }
    }
    ret = 0;
out:
    return ret;
}

static int
Simulator_init(Simulator *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int sim_ret;
    static char *kwlist[] = {"sample", "event_classes", "num_loci", 
            "num_parents", "max_population_size", "max_occupancy", 
            "dimension", "simulate_pedigree", "random_seed", "torus_diameter",
            "pixel_size", "recombination_probability", NULL}; 
    PyObject *sample, *events;
    sim_t *sim = PyMem_Malloc(sizeof(sim_t));
    self->sim = sim; 
    if (self->sim == NULL) {
        goto out;
    }
    memset(self->sim, 0, sizeof(sim_t));
    sim->num_loci = 1;
    sim->num_parents = 2;
    sim->torus_diameter = 1000;
    sim->pixel_size = 2;
    sim->recombination_probability = 0.5;
    sim->random_seed = 1;
    sim->max_population_size = 1000;
    sim->max_occupancy = 10;
    sim->dimension = 2;
    sim->simulate_pedigree = 0;
    sim->max_time = DBL_MAX;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!|IIIIIIkddd", kwlist, 
            &PyList_Type, &sample, 
            &PyList_Type, &events, 
            &sim->num_loci, &sim->num_parents, &sim->max_population_size,
            &sim->max_occupancy, &sim->dimension, &sim->simulate_pedigree,
            &sim->random_seed, &sim->torus_diameter, &sim->pixel_size, 
            &sim->recombination_probability)) {
        goto out;
    }
    if (Simulator_parse_sample(self, sample) != 0) {
        goto out;
    }
    if (Simulator_parse_events(self, events) != 0) {
        goto out;
    }
    if (Simulator_check_input(self) != 0) {
        goto out;
    }
    sim_ret = sim_alloc(self->sim);
    if (sim_ret != 0) {
        handle_library_error(sim_ret);
        goto out;
    }
    sim_ret = sim_initialise(self->sim);
    if (sim_ret != 0) {
        handle_library_error(sim_ret);
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static PyMemberDef Simulator_members[] = {
    {NULL}  /* Sentinel */
};

    

static PyObject *
Simulator_get_num_loci(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("I", self->sim->num_loci);
out:
    return ret; 
}

static PyObject *
Simulator_get_num_parents(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("I", self->sim->num_parents);
out:
    return ret; 
}

static PyObject *
Simulator_get_dimension(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("I", self->sim->dimension);
out:
    return ret; 
}

static PyObject *
Simulator_get_simulate_pedigree(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("I", self->sim->simulate_pedigree);
out:
    return ret; 
}

static PyObject *
Simulator_get_max_population_size(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("I", self->sim->max_population_size);
out:
    return ret; 
}

static PyObject *
Simulator_get_max_occupancy(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("I", self->sim->max_occupancy);
out:
    return ret; 
}

static PyObject *
Simulator_get_random_seed(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("k", self->sim->random_seed);
out:
    return ret; 
}

static PyObject *
Simulator_get_torus_diameter(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("d", self->sim->torus_diameter);
out:
    return ret; 
}

static PyObject *
Simulator_get_pixel_size(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("d", self->sim->pixel_size);
out:
    return ret; 
}

static PyObject *
Simulator_get_time(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("d", self->sim->time);
out:
    return ret; 
}


static PyObject *
Simulator_get_num_reproduction_events(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("K", 
            (unsigned long long) self->sim->num_reproduction_events);
out:
    return ret; 
}

static PyObject *
Simulator_get_recombination_probability(Simulator *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("d", self->sim->recombination_probability);
out:
    return ret; 
}

static PyObject *
Simulator_get_event_classes(Simulator *self)
{
    PyObject *ret = NULL;
    PyObject *l = NULL;
    PyObject *d = NULL;
    unsigned int j;
    event_class_t *e;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    l = PyList_New(self->sim->num_event_classes);
    if (l == NULL) {
        goto out;
    }
    for (j = 0; j < self->sim->num_event_classes; j++) {
        e = &self->sim->event_classes[j];
        d = Py_BuildValue("{s:d,s:d,s:d}", "r", e->r, "u", e->u, 
                "rate", e->rate);
        if (d == NULL) {
            goto out;
        }
        if (PyList_SetItem(l, j, d) != 0) {
            goto out;
        }
    }
    ret = l;
    l = NULL;
out:
    Py_XDECREF(l);
    return ret; 
}



static PyObject *
Simulator_individual_to_python(Simulator *self, individual_t *ind)
{
    PyObject *ret = NULL; 
    PyObject *key, *value;
    int status;
    double *x = ind->location;
    avl_node_t *node;
    int_map_value_t *imv;
    PyObject *ancestry = NULL;
    PyObject *loc = NULL;
    if (self->sim->dimension == 1) {
        loc = Py_BuildValue("d", x[0]);
    } else {
        loc = Py_BuildValue("(d,d)", x[0], x[1]);
    }
    if (loc == NULL) {
        goto out;
    }
    if (self->sim->simulate_pedigree == 1) {
        ret = loc; 
    } else {
        ancestry = PyDict_New();
        if (ancestry == NULL) {
            goto out;
        }
        for (node = ind->ancestry.head; node != NULL; node = node->next) {
            imv = (int_map_value_t *) node->item;
            key = Py_BuildValue("I", imv->key);
            if (key == NULL) {
                goto out; 
            }
            value = Py_BuildValue("I", imv->value);
            if (value == NULL) {
                Py_DECREF(key);
                goto out;
            }
            status = PyDict_SetItem(ancestry, key, value);
            Py_DECREF(key);
            Py_DECREF(value);
            if (status != 0) {
                goto out;
            }
        }
        ret = PyTuple_Pack(2, loc, ancestry);
        if (ret == NULL) {
            goto out;
        }
    }
out:
    if (self->sim->simulate_pedigree == 0) {
        Py_XDECREF(loc);
        Py_XDECREF(ancestry);
    }
    return ret; 
}
        
static PyObject *
Simulator_get_population(Simulator  *self)
{
    int err;
    unsigned int j;
    PyObject *ret = NULL;
    PyObject *l = NULL;
    PyObject *py_ind = NULL;
    avl_tree_t *pop = NULL;
    avl_node_t *node;
    uint64_t id;
    uintptr_t int_ptr;
    individual_t *ind;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    pop = PyMem_Malloc(sizeof(avl_tree_t));
    if (pop == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = sim_get_population(self->sim, pop);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    l = PyList_New(avl_count(pop));
    if (l == NULL) {
        goto out;
    }
    j = 0;
    for (node = pop->head; node != NULL; node = node->next) {
        id = *((uint64_t *) node->item);
        int_ptr = (uintptr_t) id;
        ind = (individual_t *) int_ptr;
        py_ind = Simulator_individual_to_python(self, ind);
        if (py_ind == NULL) {
            goto out;
        }
        if (PyList_SetItem(l, j, py_ind) != 0) {
            Py_DECREF(py_ind);
            goto out;
        }
        j++;
    }
    ret = l;
    l = NULL;
out:
    if (pop != NULL) {
        sim_free_population(self->sim, pop);
        PyMem_Free(pop);
    }
    Py_XDECREF(l);
    return ret; 
}


static PyObject *
Simulator_get_history(Simulator  *self)
{
    PyObject *ret = NULL;
    PyObject *pi = NULL;
    PyObject *tau = NULL;
    PyObject *pi_locus, *tau_locus;
    unsigned int j, l, n;
    int err;
    sim_t *sim = self->sim;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    if (self->sim->simulate_pedigree == 1) {
        PyErr_SetString(PyExc_NotImplementedError, 
                "Cannot get history for pedigree simulation"); 
        goto out; 
    }

    pi = PyList_New(sim->num_loci);
    if (pi == NULL) {
        goto out;
    }
    tau = PyList_New(sim->num_loci);
    if (tau == NULL) {
        goto out;
    }
    n = 2 * sim->sample_size;
    for (l = 0; l < sim->num_loci; l++) {
        pi_locus = PyList_New(n);
        if (pi_locus == NULL) {
            goto out;
        }
        err = PyList_SetItem(pi, l, pi_locus);
        if (err < 0) {
            goto out;
        }
        tau_locus = PyList_New(n);
        if (tau_locus == NULL) {
            goto out;
        }
        err = PyList_SetItem(tau, l, tau_locus);
        if (err < 0) {
            goto out;
        }
        for (j = 0; j < n; j++) {
            err = PyList_SetItem(pi_locus, j, PyLong_FromLong(sim->pi[l][j])); 
            if (err < 0) {
                goto out;
            }
            err = PyList_SetItem(tau_locus, j, 
                    PyFloat_FromDouble(sim->tau[l][j])); 
            if (err < 0) {
                goto out;
            }
        }
    }
    ret = Py_BuildValue("(O, O)", pi, tau);
out:
    Py_XDECREF(pi);
    Py_XDECREF(tau);

    return ret;
}

static PyObject *
Simulator_run(Simulator *self, PyObject *args)
{
    PyObject *ret = NULL;
    int status, not_done;
    uint64_t chunk = 8192; 
    double max_time = DBL_MAX;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "|d", &max_time)) {
        goto out;
    }
    not_done = 1; 
    self->sim->max_time = max_time;
    while (not_done) {
        status = sim_simulate(self->sim, chunk);
        if (status < 0) {
            handle_library_error(status);
            goto out;
        }
        not_done = status != 0; 
        if (PyErr_CheckSignals() < 0) {
            goto out;
        }
    }
    /* return True if complete coalescence has occured */
    ret = self->sim->time < max_time ? Py_True : Py_False;
    Py_INCREF(ret);
out:
    return ret;
}

static PyMethodDef Simulator_methods[] = {
    {"get_num_loci", (PyCFunction) Simulator_get_num_loci, METH_NOARGS, 
            "Returns the number of loci" },
    {"get_num_parents", (PyCFunction) Simulator_get_num_parents, METH_NOARGS, 
            "Returns the number of parents" },
    {"get_dimension", (PyCFunction) Simulator_get_dimension, METH_NOARGS, 
            "Returns the dimension of the simulation." },
    {"get_simulate_pedigree", (PyCFunction) Simulator_get_simulate_pedigree, METH_NOARGS, 
            "Returns 1 if we are simulating the pedigree; 0 otherwise." },
    {"get_max_population_size", (PyCFunction) Simulator_get_max_population_size, 
            METH_NOARGS, 
            "Returns the maximum size of the ancestral population"},
    {"get_max_occupancy", (PyCFunction) Simulator_get_max_occupancy, 
            METH_NOARGS, 
            "Returns the maximum occupancy of a single pixel"},
    {"get_random_seed", (PyCFunction) Simulator_get_random_seed, METH_NOARGS, 
            "Returns the random seed" },
    {"get_torus_diameter", (PyCFunction) Simulator_get_torus_diameter, 
            METH_NOARGS, "Returns the torus diameter" },
    {"get_pixel_size", (PyCFunction) Simulator_get_pixel_size, METH_NOARGS, 
            "Returns the size of a pixel" },
    {"get_time", (PyCFunction) Simulator_get_time, METH_NOARGS, 
            "Returns the current simulation time" },
    {"get_num_reproduction_events", 
            (PyCFunction) Simulator_get_num_reproduction_events, METH_NOARGS, 
            "Returns the number of reproduction events up to this point." },
    {"get_recombination_probability", 
            (PyCFunction) Simulator_get_recombination_probability, METH_NOARGS, 
            "Returns the probability of recombination between adjacent loci" },
    {"get_event_classes", (PyCFunction) Simulator_get_event_classes, METH_NOARGS, 
            "Returns the event classes" },
    {"get_population", (PyCFunction) Simulator_get_population, METH_NOARGS, 
            "Returns the state of the ancestral population" },
    {"get_history", (PyCFunction) Simulator_get_history, METH_NOARGS, 
            "Returns the history of the sample as a tuple (pi, tau)" },
    {"run", (PyCFunction) Simulator_run, METH_VARARGS, 
            "Simulates until at most the specified time. Returns True\
            if the required stopping conditions have been met and False \
            otherwise." },
    {NULL}  /* Sentinel */
};


static PyTypeObject SimulatorType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_discsim.Simulator",             /* tp_name */
    sizeof(Simulator),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)Simulator_dealloc, /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_reserved */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash  */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,        /* tp_flags */
    "Simulator objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    Simulator_methods,             /* tp_methods */
    Simulator_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Simulator_init,      /* tp_init */
};

/*===================================================================
 * IdentitySolver 
 *===================================================================
 */

static int 
IdentitySolver_check_nystrom(IdentitySolver *self) 
{
    int ret = 0;
    if (self->nystrom == NULL) {
        PyErr_SetString(PyExc_SystemError, "nystrom not initialised");
        ret = -1; 
    }
    return ret;
}



static int
IdentitySolver_parse_events(IdentitySolver *self, PyObject *py_events)
{
    int ret = -1;
    event_class_t *ec;
    int size = PyList_Size(py_events);
    if (size != 1) {
        PyErr_SetString(DiscsimInputError, "must have 1 events"); 
        goto out;
    }
    ec = PyMem_Malloc(size * sizeof(event_class_t));
    if (ec == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    ret = discsim_parse_event_classes(py_events, ec); 
    if (ret != 0) {
        goto out;
    }
    self->nystrom->event_classes = ec;
    self->nystrom->num_event_classes = size;
out:
    return ret; 
}


static int 
IdentitySolver_check_input(IdentitySolver *self) 
{
    int ret = -1;
    unsigned int j;
    nystrom_t *nystrom = self->nystrom;
    event_class_t *e;
    if (IdentitySolver_check_nystrom(self) != 0) {
        goto out;
    }
    if (nystrom->torus_diameter <= 0.0) {
        handle_input_error("must have torus_edge > 0");
        goto out;
    }
    if (nystrom->num_parents == 0) {
        handle_input_error("must have num_parents > 0");
        goto out;
    }
    if (nystrom->num_event_classes == 0) {
        handle_input_error("at least one event class required");
        goto out;
    }
    for (j = 0; j < nystrom->num_event_classes; j++) {
        e = &nystrom->event_classes[j];
        if (e->r <= 0.0 || e->r > nystrom->torus_diameter / 4.0) {
            handle_input_error("must have 0 < r < L / 4");
            goto out;
        }
        if (e->u <= 0.0 || e->u >= 1.0) {
            handle_input_error("must have 0 < u < 1");
            goto out;
        }
        if (e->rate <= 0.0) {
            handle_input_error("must have 0 < rate < 1");
            goto out;
        }
    }
    ret = 0;
out:
    return ret;
}

    

static void
IdentitySolver_dealloc(IdentitySolver* self)
{
    if (self->nystrom != NULL) {
        if (self->nystrom->event_classes != NULL) {
            PyMem_Free(self->nystrom->event_classes);
        }
        nystrom_free(self->nystrom);
        PyMem_Free(self->nystrom);
    }
    
    Py_TYPE(self)->tp_free((PyObject*)self);
}


static int
IdentitySolver_init(IdentitySolver *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int nystrom_ret;
    static char *kwlist[] = {"event_classes", "num_quadrature_points", 
            "torus_diameter", "num_parents", "mutation_rate", "max_x", 
            "integration_workspace_size", "integration_abserr", 
            "integration_relerr", NULL}; 
    PyObject *events;
    nystrom_t *nystrom = PyMem_Malloc(sizeof(nystrom_t));
    self->nystrom = nystrom; 
    if (self->nystrom == NULL) {
        goto out;
    }
    memset(self->nystrom, 0, sizeof(nystrom_t));
    nystrom->num_parents = 1;
    nystrom->torus_diameter = 100;
    nystrom->mutation_rate = 1e-6;
    nystrom->max_x = 50;
    nystrom->num_quadrature_points = 64;
    nystrom->integration_workspace_size = 100;
    nystrom->integration_abserr = 1e-6;
    nystrom->integration_relerr = 0.0;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|IdIddIdd", kwlist, 
            &PyList_Type, &events, 
            &nystrom->num_quadrature_points, &nystrom->torus_diameter, 
            &nystrom->num_parents, &nystrom->mutation_rate, 
            &nystrom->max_x, &nystrom->integration_workspace_size, 
            &nystrom->integration_abserr, &nystrom->integration_relerr)) {
        goto out;
    }
    if (IdentitySolver_parse_events(self, events) != 0) {
        goto out;
    }
    if (IdentitySolver_check_input(self) != 0) {
        goto out;
    }
    nystrom_ret = nystrom_alloc(self->nystrom);
    if (nystrom_ret != 0) {
        handle_library_error(nystrom_ret);
        goto out;
    }
    
       ret = 0;
out:
    return ret;
}

static PyMemberDef IdentitySolver_members[] = {
    {NULL}  /* Sentinel */
};

static PyObject *
IdentitySolver_solve(IdentitySolver *self)
{
    PyObject *ret = NULL;
    int nystrom_ret = 0;
    if (IdentitySolver_check_nystrom(self) != 0) {
        goto out;
    }
    nystrom_ret = nystrom_solve(self->nystrom);
    if (nystrom_ret != 0) {
        handle_library_error(nystrom_ret);
        goto out;
    }
    ret = Py_None;
    Py_INCREF(ret);
out:
    return ret; 
}

static PyObject *
IdentitySolver_interpolate(IdentitySolver *self, PyObject *args)
{
    PyObject *ret = NULL;
    double x, fx;
    if (IdentitySolver_check_nystrom(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "d", &x)) {
        goto out;
    }
    fx = nystrom_interpolate(self->nystrom, x);
    ret = Py_BuildValue("d", fx);
out:
    return ret; 
}

static PyObject *
IdentitySolver_get_torus_diameter(IdentitySolver *self)
{
    PyObject *ret = NULL;
    if (IdentitySolver_check_nystrom(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("d", self->nystrom->torus_diameter);
out:
    return ret; 
}

static PyObject *
IdentitySolver_get_integration_abserr(IdentitySolver *self)
{
    PyObject *ret = NULL;
    if (IdentitySolver_check_nystrom(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("d", self->nystrom->integration_abserr);
out:
    return ret; 
}

static PyObject *
IdentitySolver_get_integration_relerr(IdentitySolver *self)
{
    PyObject *ret = NULL;
    if (IdentitySolver_check_nystrom(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("d", self->nystrom->integration_relerr);
out:
    return ret; 
}

static PyObject *
IdentitySolver_get_max_x(IdentitySolver *self)
{
    PyObject *ret = NULL;
    if (IdentitySolver_check_nystrom(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("d", self->nystrom->max_x);
out:
    return ret; 
}


static PyObject *
IdentitySolver_get_mutation_rate(IdentitySolver *self)
{
    PyObject *ret = NULL;
    if (IdentitySolver_check_nystrom(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("d", self->nystrom->mutation_rate);
out:
    return ret; 
}


static PyObject *
IdentitySolver_get_num_parents(IdentitySolver *self)
{
    PyObject *ret = NULL;
    if (IdentitySolver_check_nystrom(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("I", self->nystrom->num_parents);
out:
    return ret; 
}

static PyObject *
IdentitySolver_get_num_quadrature_points(IdentitySolver *self)
{
    PyObject *ret = NULL;
    if (IdentitySolver_check_nystrom(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("I", self->nystrom->num_quadrature_points);
out:
    return ret; 
}

static PyObject *
IdentitySolver_get_integration_workspace_size(IdentitySolver *self)
{
    PyObject *ret = NULL;
    if (IdentitySolver_check_nystrom(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("I", self->nystrom->integration_workspace_size);
out:
    return ret; 
}



static PyMethodDef IdentitySolver_methods[] = {
    {"interpolate", (PyCFunction) IdentitySolver_interpolate, 
            METH_VARARGS, "interpolates the solution for a given x" },
    {"solve", (PyCFunction) IdentitySolver_solve, 
            METH_NOARGS, "solves the equation" },
    {"get_torus_diameter", (PyCFunction) IdentitySolver_get_torus_diameter, 
            METH_NOARGS, "Returns the torus diameter" },
    {"get_mutation_rate", (PyCFunction) IdentitySolver_get_mutation_rate, 
            METH_NOARGS, "Returns the mutation rate" },
    {"get_max_x", (PyCFunction) IdentitySolver_get_max_x, 
            METH_NOARGS, "Returns the maximum interpolation argument" },
    {"get_integration_abserr", 
            (PyCFunction) IdentitySolver_get_integration_abserr, 
            METH_NOARGS, "Returns the integration absolute error target" },
    {"get_integration_relerr", 
            (PyCFunction) IdentitySolver_get_integration_relerr, 
            METH_NOARGS, "Returns the integration relative error target" },
    {"get_num_parents", (PyCFunction) IdentitySolver_get_num_parents, 
            METH_NOARGS, "Returns the number of parents"},
    {"get_num_quadrature_points", 
            (PyCFunction) IdentitySolver_get_num_quadrature_points, 
            METH_NOARGS, "Returns the number of quadrature points"},
    {"get_integration_workspace_size", 
            (PyCFunction) IdentitySolver_get_integration_workspace_size, 
            METH_NOARGS, "Returns the GSL integration workspace size"},
   {NULL}  /* Sentinel */
};


static PyTypeObject IdentitySolverType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_discsim.IdentitySolver",             /* tp_name */
    sizeof(IdentitySolver),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)IdentitySolver_dealloc, /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_reserved */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash  */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,        /* tp_flags */
    "IdentitySolver objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    IdentitySolver_methods,             /* tp_methods */
    IdentitySolver_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)IdentitySolver_init,      /* tp_init */
};




/* Initialisation code supports Python 2.x and 3.x. The framework uses the 
 * recommended structure from http://docs.python.org/howto/cporting.html. 
 * I've ignored the point about storing state in globals, as the examples 
 * from the Python documentation still use this idiom. 
 */

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef discsimmodule = {
    PyModuleDef_HEAD_INIT,
    "_discsim",   /* name of module */
    MODULE_DOC, /* module documentation, may be NULL */
    -1,    
    NULL, NULL, NULL, NULL, NULL 
};

#define INITERROR return NULL

PyObject * 
PyInit__discsim(void)

#else
#define INITERROR return

void
init_discsim(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
    PyObject *module = PyModule_Create(&discsimmodule);
#else
    PyObject *module = Py_InitModule3("_discsim", NULL, MODULE_DOC);
#endif
    if (module == NULL) {
        INITERROR;
    }
    SimulatorType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&SimulatorType) < 0) {
        INITERROR;
    }
    Py_INCREF(&SimulatorType);
    PyModule_AddObject(module, "Simulator", (PyObject *) &SimulatorType);
    IdentitySolverType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&IdentitySolverType) < 0) {
        INITERROR;
    }
    Py_INCREF(&IdentitySolverType);
    PyModule_AddObject(module, "IdentitySolver", 
            (PyObject *) &IdentitySolverType);

    DiscsimInputError = PyErr_NewException("_discsim.InputError", NULL,
            NULL);
    Py_INCREF(DiscsimInputError);
    PyModule_AddObject(module, "InputError", DiscsimInputError);
    DiscsimLibraryError = PyErr_NewException("_discsim.LibraryError", 
            NULL, NULL);
    Py_INCREF(DiscsimLibraryError);
    PyModule_AddObject(module, "LibraryError", DiscsimLibraryError);

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}


