# Taur

Various routines and options for stripping adapters and UMIs from reads before alignment    

## Commands

```
java -jar taur.jar 
   -fastq_files "SAMPLE_ID_R1.fastq.gz;SAMPLE_ID_R2.fastq.gz" 
   -known_umi_file ./known_umis.csv 
   -umi_delim + 
   -known_umi_base_diff 1 
   -output_dir /output_path/ 
   -log_debug 
   -threads 32 
```

## Arguments

Argument | Description
---|---
fastq_files | A pair of fastq files
known_umi_file | File with list of known UMIs, header expects column 'KnownUmi'
umi_length | Used to extract a fixed number of bases from the read
umi_delim | Delimiter to use to form UMI part of read ID
known_umi_base_diff | Permitted edit distance if using a known UMI file
known_umi_use_numeric | Form UMIs using length of matched known UMI, eg READ_ID:12321:4+5
adapter_length | Adapter sequence
adapter_seq | Optionally extract an adapter sequence
output_dir | Output directory
threads | Standard for hmftools
chunk_size | 100K default, how many read groups to process per thread task
log_debug  / log_level | Standard for hmftools

