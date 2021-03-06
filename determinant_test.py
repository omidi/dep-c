import dependency_matrix 
import unittest


class KnownDeterminant(unittest.TestCase):
    def generate_matrix(self, n):    
        m = np.matrix(np.resize(np.ones(n*n, dtype=np.float64), (n,n)))
        np.fill_diagonal(m, 0)
        return m


    def test_number_of_trees(self):
        """ It test for a number of graph sizes, if the matrix determinant
        returns the correct value"""
        for n in xrange(2,15):
            m = self.generate_matrix(n)
            result = (dependency_matrix.number_of_spanning_trees(m))
            print """testing for graph size %d, and the calculated number is %d, 
            (the true value is %d) with %d difference""" % (n, round(result), n**(n-2), \
            n**(n-2) - round(result))
            self.assertEqual(float((n)**(n-2)), round(result))
            
            
    def test_number_of_forests(self):
        for n in xrange(2,15):
            m = self.generate_matrix(n)
            result = (dependency_matrix.number_of_forests(m))
            print """number of forests with %d nodes are %d""" % (n, round(result))
            # self.assertEqual(float((n)**(n-2)), round(result))


if __name__ == '__main__':
    import numpy as np
    unittest.main()
