import numpy as np

### https://doi.org/10.1002/nav.3800020109

rows = ['Alice', 'Bob', 'Carol']
cols = ['Clean', 'Sweep', 'Wash']

mat = np.array([[8, 4, 7],
                [5, 2, 3],
                [9, 4, 8]])



win = [(0,0), (1,2), (2,1)]
total = sum(mat[cell] for cell in win)
print(total)
