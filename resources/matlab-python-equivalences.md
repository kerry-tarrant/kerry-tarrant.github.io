MATLAB is presented followed by a Python equivalent. Notes end with -kt.

### ARRAYS INDICES

#### ---MATLAB---

input:
```
K=5;
NR=20;
index=0;
for j=5:NR
	for i=1:K-1
		test(i,j-4)=index;
		index=1+index;
	end
end
test
```

output:
```
test =

  Columns 1 through 9

     0     4     8    12    16    20    24    28    32
     1     5     9    13    17    21    25    29    33
     2     6    10    14    18    22    26    30    34
     3     7    11    15    19    23    27    31    35

  Columns 10 through 16

    36    40    44    48    52    56    60
    37    41    45    49    53    57    61
    38    42    46    50    54    58    62
    39    43    47    51    55    59    63
```    

#### ---PYTHON---

input:
```
import numpy as np
K=5
NR=20
ar=np.zeros((K-1,NR-4))
index=0
for j in np.arange(4,NR):
    for i in np.arange(0,K-1):
        ar[i,j-4]=index
        index=1+index
ar
print(ar)
```
output:
```
array([[ 0.,  4.,  8., 12., 16., 20., 24., 28., 32., 36., 40., 44., 48.,
        52., 56., 60.],
       [ 1.,  5.,  9., 13., 17., 21., 25., 29., 33., 37., 41., 45., 49.,
        53., 57., 61.],
       [ 2.,  6., 10., 14., 18., 22., 26., 30., 34., 38., 42., 46., 50.,
        54., 58., 62.],
       [ 3.,  7., 11., 15., 19., 23., 27., 31., 35., 39., 43., 47., 51.,
        55., 59., 63.]])
```

### MULTIPLE COMPARISON TEST

#### ---MATLAB---

C details the following

Columns 1,2: comparison

Column 3: lower CI

Column 4: stat

Column 5: upper CI

Column 6: p-val

input:
```
y = [9.87, 9.03, 6.81;
     7.18, 8.35, 7.00;
     8.39, 7.58, 7.68;
     7.45, 6.33, 9.35;
     6.41, 7.10, 9.33;
     8.00, 8.24, 8.44];
[pval,tbl,stats]=anova1(y);
C = multcompare(stats)
```
output:
```
C =

    1.0000    2.0000   -1.5318    0.1117    1.7551    0.9830
    1.0000    3.0000   -1.8618   -0.2183    1.4251    0.9367
    2.0000    3.0000   -1.9735   -0.3300    1.3135    0.8621
```

#### ---PYTHON---

C_mat details the following

Columns 0,1: comparison

Column 2: lower CI

Column 3: stat

Column 4: upper CI

Column 5: p-val

input:
```
import numpy as np
from scipy.stats import f_oneway, tukey_hsd
y = np.array([[9.87, 9.03, 6.81],
              [7.18, 8.35, 7.00],
              [8.39, 7.58, 7.68],
              [7.45, 6.33, 9.35],
              [6.41, 7.10, 9.33],
              [8.00, 8.24, 8.44]])
y0 = y[:,0].copy()
y1 = y[:,1].copy()
y2 = y[:,2].copy()
result = tukey_hsd(y0,y1,y2)
conf = result.confidence_interval(confidence_level=.95)
num_grp = len(y[0,:])
tri_num = int((num_grp-1)*(num_grp)/2)
C_mat = np.zeros((tri_num,6))
i_index = 0
j_index = 0
row=0
for i in range(num_grp):
    for j in range(num_grp):
        if i < j:
            C_mat[row,0] = i
            C_mat[row,1] = j
            C_mat[row,2] = conf.low[i,j]
            C_mat[row,3] = result.statistic[i,j]
            C_mat[row,4] = conf.high[i,j]
            C_mat[row,5] = result.pvalue[i,j]
            row+=1
print(C_mat)
```
output:
```
[[ 0.          1.         -1.5318158   0.11166667  1.75514913  0.98299484]
 [ 0.          2.         -1.8618158  -0.21833333  1.42514913  0.93674513]
 [ 1.          2.         -1.97348247 -0.33        1.31348247  0.86208746]]
```

### FLATTEN NUMPY ARRAY AND MATLAB INDEXING BY ELEMENT

#### ---MATLAB---

input:
```
yy = [1,3;
	2,4];
yy(2)
yy(4)
```
output:
```
ans = 

	2

ans = 

	4
```

#### ---PYTHON---

input:
```
import numpy as np
yy = np.array([[0,2],
		   [1,3]])
print(yy.flatten('F')[1])
print(yy.flatten('F')[3])
```
output:
```
1
3
```


### SORT, SORT INDEX
#### ---MATLAB---

input:
```
MAT = [23, 38, 11];
[S, id] = sort(MAT, 'descend')
```
output:
```
S =

    38    23    11

id = 

    2     1     3

MAT(id) = 

    38    23    11
```

#### ---PYTHON---

input:
```
import numpy as np
MAT = np.array([23, 38, 11])
S = np.sort(MAT)[::-1]
id = np.argsort(MAT)[::-1]
print('S =',S)
print('id =',id)
print('np.take(MAT, id) =',np.take(MAT, id))
```
output:
```
S = [38 23 11]
id = [1 0 2]
np.take(MAT, id) = [38 23 11]
```

### 

#### ---MATLAB---

#### ---PYTHON---
