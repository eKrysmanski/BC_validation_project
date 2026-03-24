## March 21st

Looking at using DESeq2 package as a method of normalization instead of CPM. Normalization for size should show that slope of the line is 1, and should intercept at 0,0. 
DESeq will scale the enriched and unenriched to the same relative amounts. So 1 to 1 should appear if approach is proportional. If follows this then run tost on the slope for total equilvalence across the dataset. 

## March 23rd 

Looked at the diagnostics looks pretty similar to the cpm plots. Slope of the DESeq2 normalization fit better over the 1:1 line. All of the linear models of the combined sexes and individual sexes showed equivalence to the 1:1 slope line which was expected. However it did fit TFA better than CPM normalization as that on screws with the male CPM normalization. 
Will look at creating models for individual genes then looking at if any individual genes are strangely skewed more or less so by enriching the sequencing library. 
