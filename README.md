# data-from-Peri-L-Cell-2016

Data from the paper Perié L, Cell 2016 doi: 10.1016/j.cell.2015.11.059.

MPP, HSC and CMP have been labelled with a cellular barcode. Their progenies have been sorted and their barcode analyzed. 
The data available here have been quality filtered, renormalized to 10^5 reads and average over replicates as described in the paper. The filtering code used for this is attached. It is not optimal as some steps are still done in excel but it is what I used at the time. 

The name of the text file are as follow:
"name of the experiment / type of progenitor/ time point / mouse number "

For CMP there is 3 experiments done at day 6 (LP28,31,32). The figure 1 is the average over the tree experiments. 

In the txt files, you''ll a matrix per mouse where:
- the column tag correspond to the identity of the barcodes
- the other columns correspond to the reads per barcode in a given sample. The name of the sample are coded as follow: "day_typeofprogenitor_micenumber_celltype"
M or My for myeloid cells
E for erythroblasts
Dc for dendritic cells
B for B cells

You will also find a code for plotting heatmaps in R used for this paper. 
And the script for the classification of the barcodes in categories. 
