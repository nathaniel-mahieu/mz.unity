# Tips and Troubleshooting

## Memory Limits
* The combinatorial search problem is memory intensive.  BM.limits, nrow(A), nrow(B), and nrow(M) all contribute to the increase in memory requirements.  If you have memory issues try to limit the number of peaks being searched or decrease the search depth.
* It is convenient both for memory and clarity to divide searches into groups which are related in their depth and type.  Many searches require only a depth of 1 as all the intermediate peaks are present. These can be separated from the larger depth searches such as cross-polarity and mer searches.


