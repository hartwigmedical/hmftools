select
	seq_region.name as chromosome,
    t.seq_region_start as transcript_start,
    t.seq_region_end as transcript_end,
    t.stable_id as transcript_id,
    t.version as transcript_version,
    xref.display_label as gene_name,
    gene.stable_id as gene_id,
    gene.seq_region_start as gene_start,
    gene.seq_region_end as gene_end,
    GROUP_CONCAT(karyotype.band SEPARATOR '-') as chromosome_band,
	dbprimary_acc as entrezId,
    e.stable_id as exon_id,
    e.seq_region_start as exon_start,
    e.seq_region_end as exon_end
from gene
	join object_xref on  gene.gene_id=object_xref.ensembl_id
    join xref on (xref.xref_id=object_xref.xref_id and external_db_id = 1300)
    join karyotype on gene.seq_region_id=karyotype.seq_region_id
	join transcript t on gene.canonical_transcript_id = t.transcript_id
    join seq_region on t.seq_region_id = seq_region.seq_region_id
    join exon_transcript et on et.transcript_id = t.transcript_id
    join exon e on et.exon_id = e.exon_id
where
    ((gene.seq_region_start >= karyotype.seq_region_start and gene.seq_region_start <= karyotype.seq_region_end)
		or (gene.seq_region_end >= karyotype.seq_region_start and gene.seq_region_end <= karyotype.seq_region_end))