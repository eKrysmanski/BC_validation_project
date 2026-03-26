## March 21st

Looking at using DESeq2 package as a method of normalization instead of CPM. Normalization for size should show that slope of the line is 1, and should intercept at 0,0. 
DESeq will scale the enriched and unenriched to the same relative amounts. So 1 to 1 should appear if approach is proportional. If follows this then run tost on the slope for total equilvalence across the dataset. 

## March 23rd 

Looked at the diagnostics looks pretty similar to the cpm plots. Slope of the DESeq2 normalization fit better over the 1:1 line. All of the linear models of the combined sexes and individual sexes showed equivalence to the 1:1 slope line which was expected. However it did fit TFA better than CPM normalization as that on screws with the male CPM normalization. 
Will look at creating models for individual genes then looking at if any individual genes are strangely skewed more or less so by enriching the sequencing library. 

## March 24th 
TOST individual genes, using mean and sd of normalized data, 
  Steps: 
  mean and sd of normalized values for each gene
  Perform TOST for each gene
  Extract p-values and plotting data from tests into same data frame 
  Pass fail logical test p-value <0.05 
                  nhst = significantly different from 0 
                  tost = practically equivalent within equivalence bounds 
## March 26th

Looking at the data from the individual genes and doing TOST on combined data set and sex separated data sets having no genes that passed the TOST. Will try making regressions against the linear models instead of the actual pure read data. Given that he normalization of DESeq should make the enriched and unenriched values the same the slope should equal 0. 

