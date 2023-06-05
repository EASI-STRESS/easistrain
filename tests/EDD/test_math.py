import numpy
from easistrain.EDD import math


def test_scattering_vector():
    angles = numpy.random.uniform(-180, 180, (100, 5))
    results1 = math.compute_qs(angles)
    results2 = numpy.asarray([math.qsample(a) for a in angles]).T
    numpy.testing.assert_allclose(results1, results2)
