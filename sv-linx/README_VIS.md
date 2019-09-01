# LINX Visualisation


# Parameterization
//TODO: max_distance_labels

## Relative Track Sizes

The relative sizes of the gene, segment and copy number tracks are controlled with `gene_relative_size`, `segment_relative_size` and 
`cna_relative_size` parameters respectively.

Argument | Default | Description 
---|---|---
gene_relative_size| 0.3 | Size of gene track relative to segments and copy number alterations
segment_relative_size | 1 | Size of segment track relative to copy number alterations and genes
cna_relative_size | 2 | Size of gene copy number alteration relative to genes and segments


## Font Size

The following parameters control the font size. 

Argument | Default | Description 
---|---|---
min_label_size| 35 | Minimum size of labels in pixels
max_label_size | 40 | Maximum size of labels in pixels
max_distance_labels | 100 | Maximum allowed number of distance labels
max_gene_characters | 5 | Maximum allowed gene length before applying scaling

The label size scales linearly from the min to the max label size as an inverse function of the number of distance labels to be plotted. 
If the number of distance labels exceeds `max_distance_labels`, no distance labels will be shown and all labels will be sized with `min_label_size`.

The same label size will be applied genes unless there is a gene which exceeds `max_gene_characters` in length. In this case, all genes
will be scaled down to prevent the gene labels from going outside the gene track. Adjusting this parameter is best done in conjuction with
the `gene_relative_size` parameter.


## Chromosome Range Panel

Argument | Default | Description 
---|---|---
chr_range_height| 150 | Chromosome range row height in pixels
chr_range_columns| 6 | Maximum chromosomes per row


## Fusion Panel

Argument | Default | Description 
---|---|---
fusion_height| 250 | Height of each fusion in pixels
fusion_legend_rows| 1 | Number of rows in protein domain legend
fusion_legend_height_per_row| 35 | Height of each row in protein domain legend 

