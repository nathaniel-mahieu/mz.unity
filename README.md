# mz.unity
mz.unity is an R package for detecting and exploring complex relationships in accurate-mass mass spectrometry data.  Mz-unity implements a combiunatorial search which can be tailored to many specific relationships.  These include simple relationships like isotopes, charge carriers, and common neutral losses but also complex relationships such as distal fragments, mers, and background mers.  

For more detailed information please see the forthcoming [publication](#).

## R Package
mz-unity is available as an R package on GitHub: [nathaniel-mahieu/mz-unity](https://github.com/nathaniel-mahieu/mz-unity)

## Installation
```r
#install.packages("devtools")
devtools::install_github("nathaniel-mahieu/warpgroup")
```

## Usage
```r
library(mz-unity)
relationships = mz.unity.search(A, B, M, ppm, BM.limits)
```

### Parameters
* A, B, M: matrices with accurate mass, charge, and direction (denoting gain or loss).  A is searched for in combinations of B and M.
* ppm: the maximum mass error allowed
* BM.limits: The combinatorial depth with which to search.

### Example
Toy data and more examples can be found in the [/inst directory](/inst/). Examples are copied from there.

```r
pIso2 = mz.unity.search(A = ps, 
                        B = ps, 
                        M = M.iso, ppm = 1, 
                        BM.limits = cbind(M.min = c(1), M.max = c(1), B.n = c(1)))
head(pIso2)
```

```
##      A        ppm M.1 B.1
## 1   25 -0.1904717   1   5
## 2   10  0.1236177   1  25
## 3   57  0.3110863   1  48
## 4   94 -0.8560973   1  57
## 5 1015  0.3189293   1 921
## 6 1083  0.3035952   1 938
```

Below we annotate and plot some simple relationships:

```r
rels = mz.unity.search(A = ps, 
                       B = ps, 
                       M = M.z, ppm = 10, 
                       BM.limits = cbind(M.min = c(2), M.max = c(2), B.n = c(1))
                       )
```

```r
df = subset(rels, rel %in% c("z", "pol", "nl"))
df = reindex(df, '^A$|^A.|^B$|^B.', mzz$id)

vertices = c(unique(subset(mzz, monoiso & source %in% c("psn", "psp"))[,"id"]))
vertices = unique(c(vertices, unlist(expandGraph(df)[,1:2])))

g = graph.data.frame(expandGraph(df)[,1:2], vertices = data.frame(v=as.character(vertices)))
plot(g, asp = .4, vertex.size = 2, edge.arrow.mode=1, edge.arrow.size = 0.3, vertex.label.cex=.6, vertex.frame.color='transparent', vertex.label.color = "transparent")
```

<img src="inst/figure_examples/unnamed-chunk-22-1.png" title="plot of chunk unnamed-chunk-22" alt="plot of chunk unnamed-chunk-22" style="display: block; margin: auto;" />

#### Complex relationship plot
This plot is of the relationships involving glutamate, NAD and their dimer.  Aggregation of compound groups was performed with `aggregate.self` under the assumption that all fragments were mers. This assumption is incorrect but necessary because more information is needed to distinguish fragments from mers.
<img src="inst/figure_examples/unnamed-chunk-26-2.png" title="plot of chunk unnamed-chunk-26" alt="plot of chunk unnamed-chunk-26" style="display: block; margin: auto;" />


# See Also
- [Complex Examples](/inst/examples.md)
- [Troubleshooting and Tips](/inst/troubleshooting.md)

# License
This project is licensed under the terms of the GPL-3 license.