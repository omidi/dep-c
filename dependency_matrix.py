
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
        matrix[index, index] = -np.sum(row)        
        index += 1
    return matrix
    

if __name__ == '__main__':
    import csv
    import numpy as np
    m = load_matrix('matrix_file', sep=';')
    if m.any():
        print m
        print laplacian_matrix(m)
    
