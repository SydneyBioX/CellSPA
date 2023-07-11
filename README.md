# CellSPA

Cell Segmentation Performance Assessment (CellSPA) provides an evaluation 
framework for in situ cell segmentation results.


## Installation

R package `CellSPA` can be installed as follows.

```{r}
BiocManager::install("SydneyBioX/CellSPA")
```


Example data files for XeniumBreastCancer1 dataset may be found here: https://drive.google.com/drive/folders/1hvNFDFteLp_S7qwDJgucLlLETL3f0rTs?usp=sharing. Please place under ``vignettes/data``. In this data directory, place the output .csv files from BIDCell (with format ``cell_outputs_{number}.csv``, and exclude ``cell_expr.csv``) into ``BIDCell_csv_output``. Place the segmentation .tif file (e.g, ``epoch_1_step_4000_connected_v3.tif``) in the data directory as well.

## Citation

If CellSPA has assisted you with your work, please kindly cite our paper:

Fu, X., Lin, Y., Lin, D., Mechtersheimer, D., Wang, C., Ameen, F., Ghazanfar, S., Patrick, E., Kim, J., & Yang, J. Y. H. (2023). Biologically-informed self-supervised learning for segmentation of subcellular spatial transcriptomics data. bioRxiv, 2023.2006.2013.544733. https://doi.org/10.1101/2023.06.13.544733


