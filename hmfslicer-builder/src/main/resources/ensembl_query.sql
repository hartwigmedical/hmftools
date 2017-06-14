select
	seq_region.name as chromosome,
    transcript.seq_region_start as transcript_start,
    transcript.seq_region_end as transcript_end,
    transcript.stable_id as transcript_id,
    transcript.version as transcript_version,
	xref.display_label as gene_name,
    GROUP_CONCAT(karyotype.band SEPARATOR '-') as chromosome_band,
	dbprimary_acc as entrezId
from gene
	join object_xref on  gene.gene_id=object_xref.ensembl_id
    join xref on (xref.xref_id=object_xref.xref_id and external_db_id = 1300)
    join karyotype on gene.seq_region_id=karyotype.seq_region_id
	join transcript on gene.canonical_transcript_id = transcript.transcript_id
    join seq_region on transcript.seq_region_id = seq_region.seq_region_id
where
    ((gene.seq_region_start >= karyotype.seq_region_start and gene.seq_region_start <= karyotype.seq_region_end)
		or (gene.seq_region_end >= karyotype.seq_region_start and gene.seq_region_end <= karyotype.seq_region_end))