import numpy as np

def load_matrix(infile, sep=','):
    try:        
        with open(infile) as f:
            csv_file = csv.reader(f, delimiter=sep)
            file_content = '; '.join([' '.join(map(str, row))
                                      for row in csv_file])
            try:
                m = np.matrix(file_content)
            except ValueError:
                print 'Something wrong in the content of input file!'
                return np.matrix([])
    except IOError:
        print 'IOError: perhaps file %s does not exist!' % infile
        return np.matrix([])
    return m


def laplacian_matrix(matrix):
    index = 0
    matrix *= -1
    for row in matrix:
        matrix[index, index] = -1 * np.sum(row)        
        index += 1
    return matrix


def matrix_minor(matrix, i, j):
    n = len(matrix)    
    return matrix[np.delete(np.arange(n), i),:][:, np.delete(np.arange(n), j)]


def matrix_minor_last_row_column(matrix):
    n = len(matrix) - 1
    return matrix[0:n, 0:n]


from numpy.ctypeslib import ctypes    
from numpy.ctypeslib import ndpointer    
lib = ctypes.cdll.LoadLibrary("./determinant.so")
calculate_determinant = lib.determinant
calculate_determinant.restype = ctypes.c_double
calculate_determinant.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                ctypes.c_int]

        
def number_of_spanning_trees(matrix):
    L_minor = matrix_minor_last_row_column(laplacian_matrix(matrix))
    return np.linalg.det(L_minor)
    # return np.exp(calculate_determinant(L_minor.flatten(), len(L_minor)))
    
            

if __name__ == '__main__':
    import csv
    import numpy as np
    import intdet 
    n = 9
    m = np.matrix(np.resize(np.ones(n*n, dtype=np.float64), (n,n)))
    print intdet._gauss_det(m)
    np.fill_diagonal(m, 0)
    print m
