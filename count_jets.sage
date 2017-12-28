from itertools import takewhile
import pprint
def pascal_matrix(size=10):
    '''
    Creates a matrix of size x size, whose entries are numbers in zp, such that
    entry M_{ij} is obtained recursively as in pascal triangle. Each diagonal on 
    the matrix corresponds to a row in Pascal's triangle.
    '''
    size = size
    M = ones_matrix(size,size)
    for i in range(1, size):
        for j in range(1, size):
            M[i, j] = M[i - 1, j] + M[i, j - 1]
    return M


def block_pascal_column(c, l):
    matrices = []
    M = pascal_matrix(size=c).cholesky()
    for i in range(1, l):
        R = matrix(c-i, i, 0).augment(pascal_matrix(size=c-i).cholesky())
        M = M.stack(R)
    return M

def zeros(m,n):
    return matrix(m,n,0)

def block_pascal(d, l):
    M = block_pascal_column(d,l)
    r = d
    c = d-1
    R = zeros(r, c)
    for i in range(1,l):
        A = block_pascal_column(d-i, l-i)
        M = M.augment(R.stack(A))
        r = r + c
        c = c - 1
        R = zeros(r, c)
    return M

def generate_matrices(i, n):
    l = n/3
    h=0
    ary = []
    for c in range(((n+i)+1)):
        num_rows = [0 for _ in xrange(l+1)]
        rows_to_keep = [[] for _ in xrange(l+1)]
        columns_to_keep = [[] for _ in xrange(l+1)]
        r = c+1
        t = 0
        a = 0
        for k in range(l+1):
            m = n - 3*k
            j = i + k
            d = c - k
            if j >= m + 2*k:
                continue
            if d < k-c:
                continue

            if m == 0:
                if d <= (j + m):
                    a = 1
                    r -= d
                else:
                    a = 0
            elif d in range(j+1) and d <= m:
                a = d+1
            elif d in range(j, m):
                a = j+1
            elif d in range(m, d+m):
                a = min(max(m-(d-j)+1,0), m+1)
                r = r-(d-m)
            num_rows[k] = a
            rows_to_keep[k] = range(r-a, r)
            if d in range(m, d+m) and m != 0:
                r += d-m
            if d >= 0:
                columns_to_keep[k] = range(t, t+min(d+1, max((m-j)/2,0)))
            else:
                columns_to_keep[k] = 0
            r += d
            t += d + 1

        cols= tuple([item for sublist in columns_to_keep for item in sublist])
        rows= tuple([item for sublist in rows_to_keep for item in sublist])
        pcols = list(columns_to_keep)
        prows= list(rows_to_keep)
#         for idx in range(1,len(pcols)):
#             for jdx in range(len(pcols[idx])):
#                 pcols[idx][jdx] -= ((c+1)*(c+2))/2 - ((c+1-idx)*(c+2-idx))/2
#         for idx in range(1,len(prows)):
#             for jdx in range(len(prows[idx])):
#                 prows[idx][jdx] -= ((c+1)*(c+2))/2 - ((c+1-idx)*(c+2-idx))/2
#         print rows_to_keep, columns_to_keep
        truncated_matrix = block_pascal(c+1,min(c+1,l+1))[rows,cols]
        h += truncated_matrix.rank()
        ary += [truncated_matrix.rank()]
#        print '----------------------------------'
#        print prows, pcols
#        print 'c = ', c
#        print 'rows = ', rows_to_keep
#        r_to_keep = []
#        for j in rows_to_keep:
#            if j != []:
#                r_to_keep += [len(j)]
#        print r_to_keep
#        print 'num_cols = ', [len(i) for i in takewhile(lambda x: x!= [], columns_to_keep)]
#        print 'matrix\n', truncated_matrix.str()
#        print truncated_matrix.rank()

#        print '(i,c) = (%d,%d), rank = %d' % (i,c,truncated_matrix.rank())
#        print '----------------------------------'
      
#    print h
    return h, ary
    #matrix = block_pascal(c,min(c,l))
    #print matrix.str()
sum=0

def count_jets(n):
    sum = 0
    large_ary = []
    for i in range(-n/3,n, 2):
        print 'i = ', i
        h, ary = generate_matrices(i, n)
        sum += h
        large_ary += [ary]
    return sum, large_ary
print count_jets(9)[0]
# for i in range(11, 111):
#     print 'n = (%d, %d)' % (3*i, count_jets(3*i)[0])
#block_pascal_column(5, 4).augment(matrix(5, 4, 0).stack(block_pascal_column(4, 3))

















