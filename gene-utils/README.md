# HMF Gene Utilities

This module holds code to generate static resources from ensembl. Schema of the ensembl database can be found 
on https://m.ensembl.org/info/docs/api/core/core_schema.html


### Generating cached Ensembl data files
To annotate SVs with gene information and to support fusion detection, Linx uses gene, transcript, exon and protein domain information from the Ensembl database. 
To improve performance, this data is first extracted into 4 CSV data files and then loaded into memory each time Linx runs.

To generate these 4 data files, first run Linx with these command line options:

```
java -cp gene-utils.jar com.hartwig.hmftools.geneutils.ensembl.GenerateEnsemblDataCache
    -ensembl_db [see below] -ensembl_user "anonymous" -ensembl_pass "" 
    -output_dir /path_to_write_data_files/ -ref_genome_version [37 or 38]
```

Ensembl database URLs for 19/37 & 38 are:
- mysql://ensembldb.ensembl.org:3337/homo_sapiens_core_89_37
- mysql://ensembldb.ensembl.org:3306/homo_sapiens_core_102_38

By default Linx will use version 37, but this can be overridden using the ref_genome_version config described above.

Note that ENST00000467125 is blacklisted from Ensembl as it is shares a splice boundary with a chimeric pathogenic GOPC_ROS1 fusion transcript.
