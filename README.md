# expressyouRcell
expressyouRcell is a unique and easy-to-use R package which provides an intuitive approach for visualizing and presenting multi-dimensional variations of gene expression levels across time and space. This tool gives the possibility of generating animations of pictographic representations of cells, or pictograms, providing a convenient and intuitive method for visualizing and understanding time course variations in cellular compartments. 
A range of customizable options is provided to create cellular pictograms starting from your data! :)

## 0) Before starting
### Dependencies
* CRAN
	- data.table (>= 1.13.6),
	- ggplot2 (>= 3.3.3),
	- rsvg (>= 2.1),
	- grImport2 (>= 0.2),
	- magick (>= 2.5.2),
	- ggpubr (>= 0.4.0),
	- RColorBrewer (>= 1.1-2)  
	 
* Bioconductor
	- IRanges (>= 2.24.1),
	- org.Mm.eg.db (>= 3.12.0),
	- clusterProfiler (>= 3.18.0),  
	- DOSE (>= 3.16.0),	
	- multtest (>= 2.46.0),

### Installation
You can install expressyouRcell directly from GitHub. To do so, the devtools package is required. If not already installed on your system, run:

```
install.packages("devtools")
```

Otherwise, load devtools and install expressyouRcell by:

```
library(devtools)
install_github("https://github.com/gittina/expressyouRcell", dependencies = TRUE)
```

## 1) Load the package
To load expressyouRcell run:
```
library(expressyouRcell)
```

## 2) Create the gene-localization table
This step can be skipped if you want to use your own table with information on the localization of genes within the cellular compartments and organelles. If you provide your own table, this must contain two columns: one with gene names (named exactly ```gene_symbol```) and one with the associated information of the localization of that gene within the cell (named exactly ```subcell_struct```).

Otherwise, you can create the gene-localization table with the map_gene_localization function provided within expressyouRcell. As far as now, the annotation is working only for mouse.

### map_gene_localization function
A gene annotation file, in GTF format, is required as input. On this complete set of gene symbols, a gene ontology enrichment analysis is performed to associate a gene with a term in the cellular component ontology. For this purpose, only the sub-ontology of the cellular components is taken into consideration. This function generates the gene-localization table, which maps each gene to the locations in the cellular structures, either cellular compartments or macromolecular complexes. 

```map_gene_localization``` returns a ```data.table``` containing for each gene its localization to a subcellular structure:

| gene_symbol | subcell_struct  |
| :---:       |      :-:        |
| Rgs20       | Golgi_apparatus |
| Vcpip1      | Golgi_apparatus |
| ...      | ... |
| Sox17       | nucleus         |
| Lypla1      | nucleus         |


## 3) Prepare your data
expressyouRcell is optimized for representing multiple sets of gene expression data (e.g. multiple stages). For this reason, the input has to be organized as a list of ```data.table```s. For example, in case of multiple stages, each ```data.table``` should correspond to a specific time point. 

Each ```data.table``` must have at least a column of gene names named precisely ```gene_symbol```. 
The input table can also contain additional columns with values of gene expression levels (CPM or FPKM) or results from upstream differential analysis pipeline (such as fold changes and pvalues).

| gene_symbol | Value  |
| :---:       |      :-:        |
| Rgs20       | 0.5 |
| Vcpip1      | 1.2 |
| Sox17       | -2.5         |
| Lypla1      | 1.1         |

### plot_cell function
This function simply allows the user to visualize the chosen cellular map with the default colors. The function requires as input the data.table with the graphical information (coordinates and colors for the cellular organelles). 

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
If  ```coloring_method``` is equal to ```mean``` or ```median```,  genes are first grouped according to their localization, then, mean (or median) of numeric values associated with each gene is computed for each group. In this case, the ```data.table```s  in the input list must also have an additional column containing numeric values (e.g. logFC, CPM values, etc.). To specify the value column on which you want to base your coloring, an additional parameter with the name of the column (```col_name```) must be provided as input. The given name must be compatible with column names in the ```data.table```s of the input list. 

expressyouRcell can handle your output in two main different manners, and this can be achieved with the optional parameters in the ```color_cell``` function.

#### 1) Classify genes into separate groups and for each one generate a distinct plot
expressyouRcell allows you to selectively visualize only genes belonging to distinct classes (e.g. either "up-" or "down-regulated" genes) and generating separate plots for each of the specified categories of genes. In this case, separate analysis for each subset of genes can be performed, and expressyouRcell will then output single ```ggplot``` objects for each category. 
To select this analysis, you must specify as the ```group_by``` parameter the name of a column with the categorical variable (e.g. “class”) on which you have previuoly stored the gene classification. The additional parameter ```grouping_vars``` can be specified to subselect the categories you are interested to plot (e.g. in case of DEGs classification, “up” and “down”). Default value of this parameter is null. In this case, all the genes are selected and their corresponding values are averaged for each subcellular localization, regardless any classification.

For example, the following lines will generate two distinct cellular pictograms (for the specified classes '+' and '-') for each time point in your list, as can be seen in the picture below.

```
example_list_output <- color_cell(timepoint_list = example_list,
                                  plot_data = neuron_dt_nocyto,
                                  gene_loc_table = gene_loc_table,
                                  coloring_mode = "mean",
                                  group_by = "class",
                                  grouping_vars = list("class"=c("+","-")),
                                  col_name = "logFC",
                                  colors = list("+" = c("#eaf3ea", "#307e2d"),
                                                "-" = c("#f3eaea", "#7e302d")))
```

![alt text](https://github.com/gittina/expressyouRcell/blob/master/vignettes/readme_img1.png?raw=true)

If your data have not been previously organized into distinct classes of genes, expressyouRcell can perform this step for you. To do so, you have to provide the tool with some additional parameters, such as the cutoff values for the identification of significant differentially expressed genes:
* ```thr``` to specify the cutoff value to be applied on the ```col_name``` column, 
* ```pval_col``` to specify the column containing the statistical significance values,
* ```pval_thr``` to specify the cutoff value  to be applied on the ```pval_col``` column.

#### 2) Classify genes into separate groups and merge all the results into a single plot 
If you do not want to discriminate genes by defined categories, you can set the ```group_by``` parameter to null. This is also the default value. In this case, no grouping by classification value is performed, and values of genes mapped to each subcellular localization are averaged regardless their classification and plotted together on the cellular pictograms. However, to avoid poorly informative pictograms, it is recommended to include only differentially expressed genes with the ```grouping_vars```, in particular when the logFC values are used for defining the color shade of the cellular regions.

The following lines will then output the cellular pictograms in the picture below.

```
example_list_output_together <- color_cell(timepoint_list = example_list,
                                           plot_data = neuron_dt_nocyto,
                                           gene_loc_table = gene_loc_table,
                                           coloring_mode = "mean",
                                           grouping_vars = list("class"=c("+","-")),
                                           col_name = "logFC",
                                           colors = list("+" = c("#eaf3ea", "#307e2d"),
                                                         "-" = c("#f3eaea", "#7e302d")))
```


 <img src="https://github.com/gittina/expressyouRcell/blob/master/vignettes/readme_img2.png" width="450" height="380">

### Enrichment based p-value
Enrichment analysis restricted to the sub-ontology of cellular components is performed on genes input by the user. Colors of each subcellular compartment are based on pvalues from the Fisher’s test, used to assess the statistical significance of the enrichment.

## create_animation function
Create animation between stages
