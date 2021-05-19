# expressyouRcell
expressyouRcell is a unique and easy-to-use R package which provides an intuitive approach for visualizing and presenting multi-dimensional variations of gene expression levels across time and space. This tool gives the possibility of generating animations of pictographic representations of cells, or pictograms, providing a convenient and intuitive method for visualizing and understanding time course variations in cellular compartments. 
A range of customizable options is provided to create cellular pictograms starting from your data! :)

## 1) Load the package
To load expressyouRcell run:
```
library(expressyouRcell)
```

## 2) Create the gene-localization table
This step can be skipped if you want to use your own table with information on the localization of genes within the cellular compartments and organelles. If you provide your own table, this must contain two columns: one with gene names (named exactly "gene_symbol") and one with the associated information of the localization of that gene within the cell (named exactly subcell_struct).

Otherwise, you can create the gene-localization table with the map_gene_localization function provided within expressyouRcell. As far as now, the annotation is working only for mouse.

### map_gene_localization function
A gene annotation file, in GTF format, is required as input. On this complete set of gene symbols, a gene ontology enrichment analysis is performed to associate a gene with a term in the cellular component ontology. For this purpose, only the sub-ontology of the cellular components is taken into consideration. This step generates the gene-localization table, which maps each gene to the locations in the cellular structures, either cellular compartments or macromolecular complexes. 

## 3) Prepare your data
expressyouRcell is optimized for representing multiple sets of gene expression data (e.g. multiple stages). For this reason, the input has to be organized as a list of data.tables. For example, in case of multiple stages, each data.table should correspond to a specific time point. 
Each data.table must have at least a column of gene names named precisely “gene_symbol”.

### plot_cell function
This function simply allows the user to visualize the chosen cellular structure with the default colors. The function requires as input the data.table with the graphical information (coordinates and colors for the cellular organelles). 

Example of usage:
```
plot_cell(coords_dt = neuron_dt_nocyto)
```

## color_cell
The main function is called ```color_cell``` and needs at least three mandatory input parameters.
* A list of ```data.table```s, each corresponding to a time point and must have at least a column named ```external_gene_name``` with gene symbols.
* A ```data.table``` containing ```x``` and ```y``` coordinates, name of the subcellular structures and an associated default color.
* A ```data.table``` storing for each gene the mapping to a subcellular localization according to the cellular component gene ontology. 

Different methods for assigning colors to subcellular localizations can be chosen through the ```coloring_method``` parameter. 

### Mean (or median) of values
If  ```coloring_method``` is equal to ```mean``` or ```median```,  genes are grouped according to their localization and mean (or median) of numeric values associated with each gene is computed for each group. In this case, the ```data.table```s  in the input list must also have an additional column containing numeric values (e.g. logFC, CPM values, etc.). To specify the value column on which you want to base your coloring, an additional parameter with the name of the column (```col_name```) must be provided as input. The given name must be compatible with column names in the ```data.table```s of the input list. 

expressyouRcell can handle your output in two main different manners, and this can be achieved with the optional parameters in the ```color_cell``` function.

#### 1) Classify genes into separate groups and for each one generate a distinct plot
Let's suppose you want to color only genes belonging to one class (either "up" or "down"), and obtain a separate plot for each class. In this case, separate analysis for each subset of genes can be obtained and expressyouRcell will then output single ```ggplot``` objects for each category. 
To select this analysis, you must specify as the ```group_by``` parameter the name of a column with the categorical variable (e.g. “class”) on which you have previuoly stored the gene classification. 

If you have not previously organized your genes in distinct classes, expressyouRcell can do this for you. To do so, you have to provide some additional parameters, such as the cutoff values for the identification of significant differentially expressed genes with the parameters:
* ```thr``` to specify the cutoff value on the ```col_name``` column, 
* ```pval_col``` to specify the column with the statistical significance values,
* ```pval_thr``` to specify the cutoff value on the ```pval_col``` column.

#### 2) Classify genes into separate groups and merge all the results into a single plot 
If you prefer to obtain a single plot without any discrimination of genes the ```group_by``` parameter must be null. This is also the default value. In this case, no grouping by classification is performed, and values of genes mapped to each subcellular localization are averaged regardless their classification. 

With both the options, an additional parameter ```grouping_vars``` can be specified to subselect the categories you are interested to plot (e.g. in case of DEGs classification, “up” and “down”). Default value of this parameter is null. In this case, all the genes are selected and their corresponding values are averaged for each subcellular localization, regardless any classification.

### Enrichment based p-value
Enrichment analysis restricted to the sub-ontology of cellular components is performed on genes input by the user. Colors of each subcellular compartment are based on pvalues from the Fisher’s test, used to assess the statistical significance of the enrichment.

## create_animation function
Create animation between stages
