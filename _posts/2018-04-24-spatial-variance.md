
# An Important Index underrated by GISers: Spatial Variance

by Xuezhi Cang

The spatial weighting methods and spatial cross products are widly used in the GIS and related fields,the spatial vairnce, the spatial form of variance is not well defined. This article will foucus on the defination and Python implement of the spatial variance

## 1.Defination of spatial Variance

### 1.1 Sum of Square

The consept of sum of square is from ANOVA. The meaning of the sum of square is the degree of variation of data.  <br>

Normally, the Fomula of sum of square is 
$$SS = \sum_{i=1}^n (X_i-\bar{X})^2 \tag {1}$$

### 1.2 Fomula of Variance

The purpose of variance is similar with the sum of square. The difference is that the variance is devided by the totle count of data(so we can say scaled by the totle count of data). The variance is better for measure the variation, becasue it is already scaled.<br>

Normally, the Fomula of variance is 
$$Var = {\sum_{i=1}^n (X_i-\bar{X})^2 \over N-1}\tag {2}$$

If we want to extent the variance to its spatial form, we need to figure out the spatial $\bar{X}$. Since the spatial $\bar{X}$ is hard to define,we build our spatial variance on the other form of variance.This form of variance does not need the mean value and only need data in the set.
$$Var = {\sum_{i} \sum_{i \neq j} {{(X_i-X_j})^2 \over 2} \over N(N-1)}\tag {3}$$
Does the equation ${2}$ looks familiar? Yes, this is the basic form of semi-variance, which used a lot in the kriging interpolation. The  transformation between equation ${2}$ and ${3}$ can be find in [Bachmaier and Backes 2008](https://scholar.google.com/scholar?hl=en&as_sdt=0%2C14&q=Variogram+or+semivariogram%3F+Understanding+the+variances+in+a+variogram&btnG=).

### 1.3 Spatial cross product

The definiation of spatial cross is given by [Getis 1991](https://scholar.google.com/scholar?hl=en&as_sdt=0%2C14&q=Spatial+Interaction+and+Spatial+Autocorrelation%3A+A+Cross-Product+Approach&btnG=). [Gretis 1991](https://scholar.google.com/scholar?hl=en&as_sdt=0%2C14&q=Spatial+Interaction+and+Spatial+Autocorrelation%3A+A+Cross-Product+Approach&btnG=) stated that many spatial methods are built on the same theory base, which is the spatial cross product appoach.  <br>

The basic form of spatial cross product is:
$$\Gamma = \sum_{i} \sum_{i \neq j} {W_{i,j}}{(x_i-x_j)^2}\tag{4} $$
where ${W_{i,j}}$ is the weight between $i$ and $j$.<br>

Comparing the equation $3$ and $4$, the only difference is that the spatial cross product has $W_{i,j}$ which represents the degree of impontance between $x_i$ and $x_j$. If two locations ($x_i$ and $x_j$) are close, the $W_{i,j}$ should be high. Otherwise,the $W_{i,j}$ should be low or be 0.  

For a extreme case, all the data pairs ($x_i$ and $x_j$) are treated as the 1. The spatial cross product is as the same as the sum of square So, we can see that the sum of square is the spatial cross product in a spatial case, which all the spatial weights equal to 1.



### 1.4 Spatial variance

Since we already expand the sum of square to its spatial form (spatial croee product). We can expand the variance to spatial variance.<br>

The spatial variance is:
$$Var_{spatial} = {{\sum_{i} \sum_{i \neq j} {W_{i,j}}{(x_i-x_j)^2}} \over {\sum_{i} \sum_{i \neq j} {W_{i,j}}}} \tag{5} $$ <br>

Comparing the Equation ${5}$ to Equation ${2}$, we can also find that the Equation ${2}$ is the Equation ${5}$'s spetial case, which all the weights are the same.   

# 2 Implement of spatial variance

## 2.1 data preparation

To implement the spatial variance, we first creat a spatial surface including $N*N$ (such as 900*900) points. The data is recorded in a table, which includes 3 columns: X, Y and Z. The X and Y represent the coordinators of points; the Z is the field value at this location. The  


```python
import numpy as np
import time

def gimme_mesh(n):
    minval = 0
    maxval =  n-1
    return np.meshgrid(np.linspace(minval,maxval,n), np.linspace(minval,maxval,n))

# scattered input points
N_dense = 50
x_dense,y_dense = gimme_mesh(N_dense)
x_location = np.reshape(x_dense,(N_dense*N_dense))
y_location = np.reshape(y_dense,(N_dense*N_dense))    
denpend = np.random.normal(0,10,N_dense*N_dense)# 
    
data_locations_dep_nparray = np.stack((x_location,y_location,denpend)).T
```

The data can be printed below:


```python
import pandas as pd

data_pd = pd.DataFrame(data_locations_dep_nparray)

```

There are several ways to calculate the spatial variance. Here, I compare two methods. One used a Python build-in function; the other used function in scipy. The results shows that the method using SciPy function is quicker.

First, the method used the Python build-in function is showed below:


```python
start = time.clock()

from itertools import combinations
import math
def euclidean_dist(p1,p2):
    (x1,y1),(x2,y2) = p1,p2
    return math.sqrt((x2-x1)**2+(y2-y1)**2)
def weighting_dist(data):
    x = data[:,0]
    y = data[:,1]
    points = list(zip(x,y))
    weight_list = [1/euclidean_dist(p1,p2) for p1,p2 in combinations(points,2)]
    return weight_list

def cross_product_x_part(data):
    z = data[:,2]
    points = list(zip(z,z))
    x_part = [(p1[0]-p2[0])**2 for p1,p2 in combinations(points,2)]
    return x_part

weighting  = weighting_dist(data_locations_dep_nparray)
cross_product_x = cross_product_x_part(data_locations_dep_nparray)
cross_product_whole = np.array(weighting)*np.array(cross_product_x)
cross_product_whole_sum = np.sum(cross_product_whole)
weighting_sum = np.sum(weighting)

spatial_variance = cross_product_whole_sum/weighting_sum
print ("The spatial variance is",spatial_variance )

elapsed_1 = (time.clock() - start)
print("Time used:",elapsed_1)
```

    ('The spatial variance is', 208.55579378243027)
    ('Time used:', 6.364600542318659)
    

Then, the method used the SciPy function is showed below:


```python
start = time.clock()
from scipy.spatial.distance import pdist

def weighting_dist_scipy(data):
    points_locations = data[:,0:2]
    dist = pdist(points_locations,'euclidean')
    weight_list = np.reciprocal(dist)
    return weight_list

def cross_product_x_part_scipy(data):
    Z = data[:,2]
    Y = pdist(np.reshape(Z,(-1,1)),'euclidean')
    x_part = Y**2
    return x_part

weighting_scipy  = weighting_dist_scipy(data_locations_dep_nparray)
cross_product_x_scipy = cross_product_x_part_scipy(data_locations_dep_nparray)
cross_product_whole = np.array(weighting_scipy)*np.array(cross_product_x_scipy)
cross_product_whole_sum = np.sum(cross_product_whole)
weighting_sum = np.sum(weighting_scipy)

spatial_variance = cross_product_whole_sum/weighting_sum
print ("The spatial variance is" ,spatial_variance )
elapsed_2 = (time.clock() - start )
print("Time used:" ,elapsed_2 )
```

    ('The spatial variance is', 208.55579378243027)
    ('Time used:', 0.2680237660158289)
    

# 3 Conclusion


This article introduce an important spatial index: spatial variance, which is not mentioned by a lot of spatial analysis research paper. I also impliment and compare the computation speed bwteen two methosd. The results shows that the method using scipy function is much faster than the method using build-in function.
