#
# Copyright (C) 2013 Jerome Kelleher <jerome.kelleher@ed.ac.uk>
#
# This file is part of discsim.
#
# discsim is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# discsim is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with discsim.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Coalescent simulation in continuous space under the disc replacement
model.

"""
from __future__ import division
from __future__ import print_function

import sys
import math
import random
import numbers

import ercs
import _discsim

__version__ = '1.0.0'

class Simulator(object):
    """
    Simulate the extinction/recolonisation continuum disc model
    in one or two dimensions, tracking either the genetic ancestors
    (and their genealogies at multiple loci) or the pedigree ancestors.
    If ``simulate_pedigree`` is ``True``, simulate the locations of the
    pedigree ancestors of the sample; otherwise, simulate the locations
    and ancestry of the genetic ancestors of the sample over muliple
    loci.
    
    If ``simulate_kingman`` is ``True``, a homogeneous Ancestral Recombination
    Graph is simulated on the ancestral sample after a specified number of
    generations.
    """
    def __init__(self, torus_diameter, simulate_pedigree=False, 
                 simulate_kingman=False):
        """
        Allocates a new simulator object on a torus of the specified diameter.
        """
        self.torus_diameter = torus_diameter
        self.simulate_pedigree = simulate_pedigree
        self.simulate_kingman = simulate_kingman
        self.sample = None
        self.event_classes = None
        self.num_parents = None
        self.num_loci = None
        self.recombination_probability = None
        self.arg_recombination_rate = None
        self.arg_effective_population_size = None
        self.max_occupancy = None
        self.max_population_size = None
        self.pixel_size = None
        self.dimension = None
        self.random_seed = None
        self.__simulator = None


    def __set_defaults(self):
        """
        Sets up the default values for instances that have not been
        specified.
        """
        if self.recombination_probability is None:
            self.recombination_probability = 0.5
        if self.arg_recombination_rate is None:
            self.arg_recombination_rate = 1e-8
        if self.arg_effective_population_size is None:
            self.arg_effective_population_size = 20000
        if self.num_loci is None:
            self.num_loci = 1
        if self.random_seed is None:
            self.random_seed = random.randint(0, 2**31)
        if self.max_population_size is None:
            self.max_population_size = 1000
        if self.pixel_size is None:
            self.pixel_size = 2.0
        if self.num_parents is None:
            if self.simulate_pedigree:
                self.num_parents = 2
            else:
                self.num_parents = 1 if self.num_loci is 1 else 2
        if self.max_occupancy is None:
            r = self.event_classes[0].r
            u = self.event_classes[0].u
            s = self.pixel_size
            N = self.num_parents / u
            area = s + 2 * r
            if self.dimension == 2:
                area = s**2 + 4 * r * s + math.pi * r**2
            # a density of N per unit area is probably overkill
            self.max_occupancy = int(area * N)

    def __convert_sample(self):
        """
        Coverts the sample of locations in the sample instance variable to
        the format used by _ercs, a zero-indexed list. The zero'th element
        of this list must be None
        """
        if self.sample is None:
            raise ValueError("sample must be specified")
        if len(self.sample) < 2:
            raise ValueError("At least one sample must be specified")
        if self.sample[0] is not None:
            raise ValueError("zeroth element of list samples must be None")
        sample = self.sample[1:]
        last_dim = None
        dim = None
        for x in sample:
            if isinstance(x, numbers.Number):
                dim = 1
            else:
                dim = 2
            if last_dim is None:
                last_dim = dim
            if dim != last_dim:
                raise ValueError("Sample must be either 1 or 2 dimensional")
        self.dimension = dim
        return sample

    def __convert_events(self):
        """
        Converts the events to the required low-level representation.
        """
        if self.event_classes is None:
            raise ValueError("event_classes must be specified")
        if len(self.event_classes) == 0:
            raise ValueError("At least one event class must be specified")
        for ec in self.event_classes:
            if not isinstance(ec, ercs.DiscEventClass):
                raise ValueError("Only events from the disc model are supported")
        llec = [ec.get_low_level_representation() for ec in self.event_classes]
        return llec


    def __allocate_simulator(self):
        """
        Allocates a new simulator instance with the required instance
        variables.
        """
        sample = self.__convert_sample()
        events = self.__convert_events()
        self.__set_defaults()
        self.__simulator = _discsim.Simulator(sample, events,
                    num_loci=self.num_loci, torus_diameter=self.torus_diameter,
                    pixel_size=self.pixel_size, random_seed=self.random_seed,
                    recombination_probability=self.recombination_probability,
                    arg_recombination_rate = self.arg_recombination_rate, 
                    arg_effective_population_size = self.arg_effective_population_size,
                    num_parents=self.num_parents,
                    simulate_pedigree=int(self.simulate_pedigree),
                    simulate_kingman=int(self.simulate_kingman),
                    max_population_size=self.max_population_size,
                    max_occupancy=self.max_occupancy, dimension=self.dimension)

    def get_ll_object(self):
        """
        Returns the low-level simulator object.
        """
        return self.__simulator

    def run(self, until=None):
        """
        Runs the simulation until coalescence or the specified time is
        exceeded. If ``until`` is not specified simulate until complete
        coalescence. Returns True if the sample coalesced, and False
        otherwise.

        :param until: the time to simulate to.
        :type until: float
        :return: ``True`` if the sample has completely coalesced; ``False``
            otherwise.
        :rtype: Boolean
        :raises: :exc:`_discsim.LibraryError` when the C library encounters an
            error
        """
        if self.__simulator == None:
            self.__allocate_simulator()
        dbl_max = sys.float_info.max
        t = dbl_max if until == None else until
        return self.__simulator.run(t)

    def get_time(self):
        """
        Returns the current time of the simulator.
        """
        return self.__simulator.get_time()

    def get_num_reproduction_events(self):
        """
        Returns the number of reproduction events since the beginning
        of the simulation.
        """
        return self.__simulator.get_num_reproduction_events()

    def get_population(self):
        """
        Returns the current state of the population. For a pedigree simulation,
        this returns the current locations of all individuals; in a genetic
        simulation, this also returns the ancestral material mappings for
        each individual.

        :return: the current state of the population.
        :rtype: A list describing the state of each extant ancestor. For a
            pedigree simulation, this is a list of locations. For a genetic
            simulation, this is a list of tuples ``(x, a)``, where ``x``
            is the location of the ancestor and ``a`` is its ancestry. The
            ancestry is a dictionary mapping a locus to the node it occupies
            in the genealogy for that locus.
        """
        return self.__simulator.get_population()

    def get_history(self):
        """
        Returns the history of the current ancestral population. This is
        not defined for a pedigree simulation.

        :return: the simulated history of the sample, (pi, tau)
        :rtype: a tuple ``(pi, tau)``; ``pi`` is a list of lists of integers,
            and ``tau`` is a list of lists of doubles
        :raises: :exc:`NotImplementedError` if called for a pedigree
            simulation.
        """
        return self.__simulator.get_history()

    def reset(self):
        """
        Resets the simulation so that we can perform more replicates. This must
        be called if any attributes of the simulation are changed; otherwise,
        these changes will have no effect.
        """
        self.__simulator = None


    def print_state(self):
        """
        Prints a short summary of the state of the simulator.
        """
        print("torus_diameter = ", self.torus_diameter)
        print("simulate_pedigree = ", self.simulate_pedigree)
        print("simulate_kingman = ", self.simulate_kingman)
        print("sample = ", len(self.sample) - 1)
        print("event_classes = ", len(self.event_classes))
        for ec in self.event_classes:
            print("\t", ec.get_low_level_representation())
        print("num_parents = ", self.num_parents)
        print("num_loci = ", self.num_loci)
        print("recombination_probability = ", self.recombination_probability)
        print("arg_effective_population_size = ", self.arg_effective_population_size)
        print("arg_recombination_rate = ", self.arg_recombination_rate)
        print("max_occupancy = ", self.max_occupancy)
        print("max_population_size = ", self.max_population_size)
        print("pixel_size = ", self.pixel_size)
        print("dimension = ", self.dimension)
        print("random_seed = ", self.random_seed)
        n = 0 if self.__simulator is None else len(self.get_population())
        print("population size = ", n)


