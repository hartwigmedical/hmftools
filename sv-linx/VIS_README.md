# LINX Visualisation


# Parameterization


## Relative Track Sizes

The relative sizes of the gene, segment and copy number tracks are controlled with `gene_relative_size`, `segment_relative_size` and 
`cna_relative_size` parameters respectively.


## Font Size

The following parameters control the font size. 

Argument | Default | Description 
---|---|---
min_label_size| 35 | Minimum size of labels in pixels
max_label_size | 40 | Maximum size of labels in pixels
max_distance_labels | 100 | Maximum allowed number of distance labels

The label size scales linearly from the min to the max label size as an inverse function of the number of distance labels to be plotted. 
If the number of distance labels exceeds `max_distance_labels`, no distance labels will be shown and all labels will be sized with `min_label_size`.



