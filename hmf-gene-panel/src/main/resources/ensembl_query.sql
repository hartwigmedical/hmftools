select
	seq_region.name as chromosome,
    gene.seq_region_start as gene_start,
    gene.seq_region_end as gene_end,
    gene.stable_id as gene_id,
    display_xref.display_label as gene_name,
	entrez_xref.dbprimary_acc as entrezId,
    GROUP_CONCAT(karyotype.band ORDER BY karyotype.band SEPARATOR '-') as chromosome_band,
    t.stable_id as transcript_id,
    t.version as transcript_version,
    t.seq_region_start as transcript_start,
    t.seq_region_end as transcript_end,
    e.stable_id as exon_id,
    e.seq_region_start as exon_start,
    e.seq_region_end as exon_end
from gene
	join object_xref on gene.gene_id=object_xref.ensembl_id and object_xref.ensembl_object_type = 'GENE'
    join xref as entrez_xref on (entrez_xref.xref_id=object_xref.xref_id and entrez_xref.external_db_id = 1300)
    join xref as display_xref on (display_xref.xref_id=gene.display_xref_id)
    join karyotype on gene.seq_region_id=karyotype.seq_region_id
	join transcript t on gene.canonical_transcript_id = t.transcript_id
    join seq_region on t.seq_region_id = seq_region.seq_region_id
    join exon_transcript et on et.transcript_id = t.transcript_id
    join exon e on et.exon_id = e.exon_id
where
    seq_region.name in ('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT') and
    ((gene.seq_region_start >= karyotype.seq_region_start and gene.seq_region_start <= karyotype.seq_region_end)
		or (gene.seq_region_end >= karyotype.seq_region_start and gene.seq_region_end <= karyotype.seq_region_end))