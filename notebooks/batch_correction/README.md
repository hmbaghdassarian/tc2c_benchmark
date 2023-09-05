**The purpose of these analyses is to determine how the tensor decomposition used by Tensor-cell2cell is affected by batch effects between samples.** 

For use with Tensor-cell2cell, we want a dataset that represents >2 contexts. This means that we will be analyzing multiple samples, each of which have variation in the data that are due to technical factors. [Batch correction](https://doi.org/10.1038/s41592-018-0254-1) removes technical variation while preserving biological variation between samples. Batch correction is an important consideration since Tensor-cell2cell considers multiple samples to extract context-dependent patterns, and we want to make sure we are capturing true biological signals rather than sample-specific differences due to technical variability. 

Of note, datasets that contain [replicates](https://doi.org/10.1038/nmeth.3091) can be particularly useful on shedding light on batch effects when using Tensor-cell2cell. Replicates will allow us to ensure that the output factors are not simply due to technical effects. We can reasonably assume that the biological variation in samples between different contexts will be greater than that of those within contexts after using appropriate batch correction to remove technical variation. Thus, we expect Tensor-cell2cell to capture overall communication trends differing between contexts and can then assess that output factors aren't simply due to technical effects by checking that the sample loading values in the factors' context dimension are similar for biological replicates and do not have  high loadings for just one sample in the context dimension. 

However, it is not always possible to have replicates, particularly technical replicates. Thus, it is important to understand the extent to which batch effects may change Tensor-cell2cell's results and how to appropriately handle them. Additional considerations:

* 1\.  CCC inference methods typically require a gene by cell counts matrix. Many batch correction methods only return a latent space, so it is important to select a batch correction method that can return a corrected counts matrix (see Table 1 [here](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1850-9/tables/1) or [here](https://doi.org/10.1093/nargab/lqac022)).  
* 2\. Non-negativity: We typically want a non-negative counts matrix as input both for CCC inference (2a) and Tensor-cell2cell (2b). Most batch correction methods that do return a corrected counts matrix have negative values. Thus, if using one of these methods, one must typically replace negative values with 0. <br>
    * 2a\. Most CCC tools require non-negative counts as input for their scoring functions. Take the simple case of the product between a ligand and receptor. A negative ligand expression value and a negative receptor count value indicates both are lowly expressed, however their product would be a positive value. If retaining negative counts in the expression matrix, we recommend using a CCC tool that implements a scoring function which can handle this issue (e.g., an additive function). 
    * 2b\. Even if the CCC tool can appropriately handle negative counts, Tensor-cell2cell requires non-negative values for decomposition. Thus, if output communication scores are negative, Tensor-cell2cell will fill these negative counts with 0 or NaN (see the [Missing Indices](../missing_indices/README.md) discussion for more details on these fill values). By default, if the tensor contains negative values, Tensor-cell2cell will treat these as 0 during decomposition. 

To this end, we wanted to assess the need and utility of batch correction on Tensor-cell2cell's outputs. We chose two batch correction methods, [scVI](https://doi.org/10.1038/s41592-018-0229-2) and [scanorama](https://doi.org/10.1038/s41587-019-0113-3). Both methods can output a corrected counts matrix and have been [benchmarked](https://doi.org/10.1038/s41592-021-01336-8) to work well. scVI uses a machine-learning approach and outputs non-negative values, whereas Scanorama uses a MNN approach and does not output non-negative values. 

For this benchmarking, we asked the following questions:
* A) At what level of batch severity are results affected? We assess batch severity using two [metrics](https://doi.org/10.1038/s41587-020-00748-9), a "clusterability -- conservation of biological variance by cell label" based metric (NMI) and a "mixability -- removal of batch effects" based metrick ([kBET](https://doi.org/10.1038/s41592-018-0254-1)).
* B) When batch severity does affect results, does batch correction mitigate these issues?ß
* C) To what extent does the introduction of negative counts effect the output? At some point in the CCC pipeline, these values will yield 0s (either by being replaced in the counts matrix or in the tensor). It may be reasonable to assume that since negative counts represent lower expression, replacing them with 0 should not drastically affect results. On the other hand, if 0 in the corrected counts matrix does not represent a biologically low value, perhaps replacing negative counts with 0 distorts results. 