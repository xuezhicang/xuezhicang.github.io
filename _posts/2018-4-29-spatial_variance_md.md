
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
