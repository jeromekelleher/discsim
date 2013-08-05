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
    
    def test_random_values(self):
        num_tests = 1000
        for j in range(num_tests):
            pixel_size = random.uniform(1, 3)
            torus_diameter = 1024 * pixel_size
            rho = random.uniform(0, 1)
            random_seed = random.randint(0, 2**32)
            num_parents = random.randint(1, 5)
            num_loci = random.randint(0, 100)
            sample_size = random.randint(2, 50)
            max_occupancy = random.randint(sample_size, 100)
            max_population_size = random.randint(sample_size, 100)
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
                    max_occupancy=max_occupancy, random_seed=random_seed)
            self.assertEqual(s.get_num_parents(), num_parents)
            self.assertEqual(s.get_num_loci(), num_loci)
            self.assertEqual(s.get_max_population_size(), max_population_size)
            self.assertEqual(s.get_max_occupancy(), max_occupancy)
            self.assertEqual(s.get_random_seed(), random_seed)
            self.assertEqual(s.get_torus_diameter(), torus_diameter)
            self.assertEqual(s.get_pixel_size(), pixel_size)
            self.assertEqual(s.get_recombination_probability(), rho)


            

            

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

