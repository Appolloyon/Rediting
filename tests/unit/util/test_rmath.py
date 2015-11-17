from unittest import TestCase
from unittest import skip

from rediting.util import rmath

class TestMath(TestCase):

    def test_calc_mean_100(self):
        seq = '1111111111'
        self.assertEqual(rmath.calc_mean(seq), 1.0)

    def test_calc_mean_50(self):
        seq = '1111100000'
        self.assertEqual(rmath.calc_mean(seq), 0.5)

    def test_calc_mean_0(self):
        seq = '0000000000'
        self.assertEqual(rmath.calc_mean(seq), 0.0)

    def test_pearson(self):
        xvals = [43,21,25,42,57,59]
        yvals = [99,65,79,75,87,81]
        xmean = rmath.calc_mean(xvals)
        ymean = rmath.calc_mean(yvals)
        PC = rmath.calc_pearson(xvals,yvals,xmean,ymean)
        self.assertAlmostEqual(PC,0.5298,places=4)

    def test_calc_tvalue(self):
        tval = rmath.calc_tvalue(0.5298,6)
        self.assertAlmostEqual(tval,1.249,places=3)

    @skip("skipped for now")
    def test_ispolyTpercent(self):
        pass
