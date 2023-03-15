**The purpose of these analyses is to determine how the tensor decomposition used by Tensor-cell2cell handles missing values.** 

Cell-cell communication (CCC) inference tools output the following for each sample: communication scores associated with sender and receiver cell type pairs and ligand-receptor (LR) pairs (see this [Review](https://doi.org/10.1038/s41576-020-00292-x) for details). 

Tensor-cell2cell structures the output of these CCC tools such that many samples are "concatenated" into a multi-dimensional array ("tensor construction"). Consequently, cell types or LR pairs that are not present across all samples will result in missing values (at the "sample - sender cell type - receiver cell type - LR pair" tensor coordinate or index) in the samples they were absent from. Sample-specific "missing" data may be the result of any of the following:

1)  Technical limitations in measuring a gene or cell. 
2)  Computational pipelines (data processing, negative expression counts or communication scores, thresholding parameters of CCC tools, etc.) that result in the exclusion of certain measurements.
3) The cell type or LR pair is truly absent from that sample due to biological reasons (a "true biological zero"). 

Tensor-cell2cell's decomposition will handle these missing values differently depending on how they are filled in during tensor construction. Within the tool, this is handled by the `how`, `cell_fill`, `lr_fill`, and `outer_fraction` parameters. See the protocol's [Tutorial 03](INSERTLINKHERE) for details. 

If the missing values are filled with NaN (the default option in Tensor-cell2cell), Tensor-cell2cell's non-negative canonical polyadic decomposition will "[mask](http://tensorly.org/stable/user_guide/sparse_backend.html)" these values during decomposition. During iteration of the Alternating Least Squares algorithm, masked indices are randomly initialized then updated in each iteration. Specifically, the masked indices in the full tensor are updated with those imputed from the previous iteration, leading to a new optimization problem and new output set of masked values in each iteration. Conceptually, this is essentially an imputation of missing values. If missing values are filled with a floating point value, they will not be masked and be considered as such in the decomposition. For example, filling missing indices with 0 will cause these indices to be treated as "true biological zeroes". 

Thus, our goal in these analyses is to determine the effect that the fraction of missing indices, as well as the the value with which they are filled, has on decomposition results. To do so, we simulate CCC across multiple samples and construct a gold-standard tensor with no missing indices. Next, we randomly generate missing indices in the tensor and fill them either with NaN or 0. Finally, we compare the similarity of decomposition outputs between the tensor with missing indices and the gold-standard using the [CorrIndex](https://doi.org/10.1016/j.sigpro.2022.108457) metric. 

Our expectations are as follows: 

1) If Tensor-cell2cell is appropriately robust to missing indices due to technical limiations or computational pipelines that exclude measurements, similarity between the gold-standard tensor decomposition output and that of the tensor with missing indices should be high even at a high fraction of missing indices (i.e., when filled with NaN, Tensor-cell2cell should be able to accurately impute the data).
2) If Tensor-cell2cell is appropriately sensitive to missing indices because those are truly absent in the sample, similarity between the gold-standard tensor decomposition output and that of the tensor with missing indices should be low as the fraction of missing indices increases (i.e., when filled with 0, Tensor-cell2cell should be able to accurately distinguish between the gold-standard tensor and that with true biological zeroes).

These analyses, as well as considerations of what the cause of missing indices is, should guide users in their choice of associated parameter values for Tensor-cell2cell. 