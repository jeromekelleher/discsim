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
Test cases for the discsim module.
"""
from __future__ import division
from __future__ import print_function 

import unittest
import random
import optparse

import _discsim


class TestInitialiser(unittest.TestCase):
    """
    Test the initialisation code for the low level interface
    """
    
    def test_bad_events(self):
        sample = [(0,0), (0,0)]
        def f(x):
            _discsim.Simulator(sample, x, torus_diameter=10)
        self.assertRaises(_discsim.InputError, f, []) 
        self.assertRaises(_discsim.InputError, f, [{}])
        self.assertRaises(_discsim.InputError, f, [[]])
        self.assertRaises(_discsim.InputError, f, [None])
        e = {"r":None, "u":None, "rate":None}
        self.assertRaises(_discsim.InputError, f, [e])
        e = {"r":None, "u":0.5, "rate":1.0}
        radii = ["12.2", 0, -1, 10, 2.6]
        for r in radii:
            e["r"] = r 
            self.assertRaises(_discsim.InputError, f, [e])
        
        impacts = ["12", -1, 0, 1, 100]
        e["r"] = 1
        for u in impacts: 
            e["u"] = u
            self.assertRaises(_discsim.InputError, f, [e])
        e["u"] = 0.5
        for rate in ["a", None, -1000, 0]:
            e["rate"] = rate
            self.assertRaises(_discsim.InputError, f, [e])

    def test_bad_parameters(self):
        sample = [(0,0), (0,0)]
        events = [{"r":0.5, "u":0.5, "rate":0.5}]
        def f(**kwargs):
            _discsim.Simulator(sample, events, **kwargs)
        self.assertRaises(_discsim.InputError, f, torus_diameter=-1)
        self.assertRaises(_discsim.InputError, f, torus_diameter=0)
        self.assertRaises(_discsim.InputError, f, num_loci=0)
        self.assertRaises(_discsim.InputError, f, num_parents=0)
        self.assertRaises(_discsim.InputError, f, max_population_size=0)
        self.assertRaises(_discsim.InputError, f, max_occupancy=0)
        self.assertRaises(_discsim.InputError, f, max_time=0)
        self.assertRaises(_discsim.InputError, f, max_time=-1)
        self.assertRaises(_discsim.InputError, f, 
                recombination_probability=-1)
        self.assertRaises(_discsim.InputError, f, 
                recombination_probability=1.1)
        self.assertRaises(_discsim.InputError, f, 
                recombination_probability=10)
        self.assertRaises(_discsim.InputError, f, pixel_size=-1)
        self.assertRaises(_discsim.InputError, f, pixel_size=0)
        self.assertRaises(_discsim.InputError, f, pixel_size=1.1, 
                torus_diameter=10)
        self.assertRaises(_discsim.InputError, f, pixel_size=0.33, 
                torus_diameter=1)
        self.assertRaises(_discsim.InputError, f, pixel_size=2, 
                torus_diameter=1)
        self.assertRaises(_discsim.InputError, f, pixel_size=0.5, 
                torus_diameter=1)
        self.assertRaises(_discsim.InputError, f, pixel_size=1, 
                torus_diameter=3)
    
    def test_bad_sample(self):
        events = [{"r":0.5, "u":0.5, "rate":0.5}]
        def f(sample, **kwargs):
            _discsim.Simulator(sample, events, **kwargs)
        self.assertRaises(_discsim.InputError, f, [])
        self.assertRaises(_discsim.InputError, f, [""])
        self.assertRaises(_discsim.InputError, f, [()])
        self.assertRaises(_discsim.InputError, f, [(1,)])
        self.assertRaises(_discsim.InputError, f, [(10,10)], torus_diameter=1)
        self.assertRaises(_discsim.InputError, f, [(-1,-1)], torus_diameter=1)
        self.assertRaises(_discsim.InputError, f, [(-1,0)], torus_diameter=1)
        self.assertRaises(_discsim.InputError, f, [(0,1)], torus_diameter=1)

    def test_random_values(self):
        num_tests = 10 
        for j in range(num_tests):
            pixel_size = random.uniform(1, 3)
            torus_diameter = 64 * pixel_size
            rho = random.uniform(0, 1)
            random_seed = random.randint(0, 2**32)
            num_parents = random.randint(1, 5)
            num_loci = random.randint(1, 100)
            sample_size = random.randint(2, 50)
            max_occupancy = random.randint(sample_size, 100)
            max_population_size = random.randint(2 * sample_size, 100)
            max_time = random.uniform(1, 1e6)
            sample = []
            for k in range(sample_size):
                x = random.uniform(0, torus_diameter)
                y = random.uniform(0, torus_diameter)
                sample.append((x, y))
            events = [{"r":random.uniform(0.1, 10), 
                "u":random.uniform(0, 1), 
                "rate":random.uniform(0, 1000)}]
            s = _discsim.Simulator(sample, events, num_loci=num_loci, 
                    torus_diameter=torus_diameter, pixel_size=pixel_size,
                    recombination_probability=rho, num_parents=num_parents,
                    max_population_size=max_population_size, 
                    max_occupancy=max_occupancy, random_seed=random_seed, 
                    max_time=max_time)
            self.assertEqual(s.get_num_parents(), num_parents)
            self.assertEqual(s.get_num_loci(), num_loci)
            self.assertEqual(s.get_max_population_size(), max_population_size)
            self.assertEqual(s.get_max_occupancy(), max_occupancy)
            self.assertEqual(s.get_random_seed(), random_seed)
            self.assertEqual(s.get_torus_diameter(), torus_diameter)
            self.assertEqual(s.get_pixel_size(), pixel_size)
            self.assertEqual(s.get_max_time(), max_time)
            self.assertEqual(s.get_recombination_probability(), rho)
            pop = s.get_population()
            locations = [x for x, a in pop]
            ancestry = [a for x, a in pop]
            self.assertEqual(len(pop), sample_size)
            for x in sample:
                self.assertTrue(x in locations)
            for a in ancestry:
                self.assertEqual(len(a), num_loci)
                d = {}
                for v in a.values():
                    if v not in d:
                        d[v] = 0
                    d[v] += 1
                self.assertEqual(len(d), 1)
            pi, tau = s.get_history() 
            self.assertEqual(num_loci, len(pi))
            self.assertEqual(num_loci, len(tau))
            for l in range(num_loci):
                self.assertEqual(2 * sample_size, len(pi[l]))
                self.assertEqual(2 * sample_size, len(tau[l]))
                for a, t in zip(pi[l], tau[l]):
                    self.assertEqual(a, 0)
                    self.assertEqual(a, 0.0)
                
            self.assertEqual(s.get_time(), 0.0)
            self.assertEqual(s.get_num_reproduction_events(), 0)

class TestSimulation(unittest.TestCase):
    """
    Tests the simulation method to see if it works correctly.
    """
    
    def setUp(self):
        
        events = [{"r":1, "u":0.99, "rate":1}]
        L = 40
        self._sample = [(0,0), (L/2, L/2)]
        self._simulator = _discsim.Simulator(self._sample, events, 
                torus_diameter=L, num_loci=1) 
        
    def test_one_event(self):
        """
        Test if the effect of one event is correct 
        """
        status = self._simulator.run(1)
        self.assertFalse(status)
        self.assertTrue(self._simulator.get_time() > 0.0)
        self.assertEqual(self._simulator.get_num_reproduction_events(), 1)
        pop = self._simulator.get_population()
        locations = [x for x, a in pop]
        # One of the sample locations must have stayed the same and 
        # the other changed.
        s = set(self._sample + locations)
        self.assertEqual(len(s), 3)
        ancestry = [a for x, a in pop]
        self.assertEqual(len(pop), 2)
        for a in ancestry:
            self.assertEqual(len(a), 1)
            d = {}
            for v in a.values():
                if v not in d:
                    d[v] = 0
                d[v] += 1
            self.assertEqual(len(d), 1)
   
    def test_num_events(self):
        """
        Tests to ensure that the number of events run is exactly what 
        is specified.
        """
        s = 0
        for j in range(1, 6):
            n = 10**j
            status = self._simulator.run(n)
            s += n
            if not status:
                self.assertEqual(self._simulator.get_num_reproduction_events(), s)
           

    def check_single_locus_history(self, pi, tau):
        self.assertEqual(len(pi), 1)
        self.assertEqual(len(tau), 1)
        self.assertEqual(len(pi[0]), 4)
        self.assertEqual(len(tau[0]), 4)
        self.assertTrue(tau[0][-1] > 0.0)
        self.assertTrue(all(a == 0 for a in tau[0][:3]))
        self.assertEqual(pi[0][-1], 0)
        self.assertEqual(pi[0][0], 0)
        self.assertEqual(pi[0][1], 3)
        self.assertEqual(pi[0][2], 3)
         

    def test_coalescence(self):
        """
        Run's the simulation until coalescence.
        """
        status = self._simulator.run()
        self.assertTrue(status)
        pi, tau = self._simulator.get_history()
        self.check_single_locus_history(pi, tau) 

    def test_pixel_size(self):
        L = 10
        sample = [(0,0), (L/2, L/2)]
        for r in [0.5, 1, 2.33, 2.5]:
            for s in [0.25, 0.5, 1, 2, 2.5]:
                events = [{"r":r, "u":0.99, "rate":1}]
                sim = _discsim.Simulator(sample, events, 
                        torus_diameter=L, pixel_size=s) 
                self.assertTrue(sim.run()) 
                pi, tau = sim.get_history()
                self.check_single_locus_history(pi, tau) 



if __name__ == "__main__":
    usage = "usage: %prog [options] "
    parser = optparse.OptionParser(usage=usage) 
    parser.add_option("-s", "--random-seed", dest="random_seed",
            help="Random seed", default=1)
    parser.add_option("-n", "--name-case", dest="name",
            help="Run this specified test", default="test")
    parser.add_option("-i", "--iterations", dest="iterations",
            help="Repeat for i iterations", default="1")
    (options, args) = parser.parse_args()
    iterations = int(options.iterations)
    random.seed(int(options.random_seed))
    testloader = unittest.TestLoader()
    suite = testloader.loadTestsFromName(options.name)
    for i in range(iterations):
        unittest.TextTestRunner(verbosity=2).run(suite)

