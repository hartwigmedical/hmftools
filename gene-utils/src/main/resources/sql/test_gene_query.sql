# Test Ensemble query

(select
	seq_region.name as chromosome,
    gene.seq_region_start as gene_start,
    gene.seq_region_end as gene_end,
    gene.stable_id as gene_id,
    display_xref.display_label as gene_name,
	GROUP_CONCAT(DISTINCT entrez_xref.dbprimary_acc ORDER BY entrez_xref.dbprimary_acc SEPARATOR ',') as entrezId,
    GROUP_CONCAT(DISTINCT karyotype.band ORDER BY karyotype.band SEPARATOR '-') as chromosome_band,
    t.stable_id as transcript_id,
    t.version as transcript_version,
    t.seq_region_start as transcript_start,
    t.seq_region_end as transcript_end,
    e.stable_id as exon_id,
    e.seq_region_start as exon_start,
    e.seq_region_end as exon_end,
    t.seq_region_strand as strand,
	if(t.seq_region_strand = -1, ce.seq_region_end - tl.seq_end + 1, cs.seq_region_start + tl.seq_start - 1) as coding_start,
	if(t.seq_region_strand = -1, cs.seq_region_end - tl.seq_start + 1, ce.seq_region_start + tl.seq_end - 1) as coding_end
from gene
	inner join object_xref on gene.gene_id=object_xref.ensembl_id and object_xref.ensembl_object_type = 'GENE'
    inner join xref as display_xref on (display_xref.xref_id=gene.display_xref_id)
    inner join karyotype on gene.seq_region_id=karyotype.seq_region_id
	inner join transcript t on gene.canonical_transcript_id = t.transcript_id
    inner join seq_region on t.seq_region_id = seq_region.seq_region_id
    inner join exon_transcript et on et.transcript_id = t.transcript_id
    inner join exon e on et.exon_id = e.exon_id
    left join xref as entrez_xref on (entrez_xref.xref_id=object_xref.xref_id and entrez_xref.external_db_id = 1300)
    left join translation tl on tl.transcript_id = t.transcript_id
    left join exon cs on cs.exon_id = tl.start_exon_id
    left join exon ce on ce.exon_id = tl.end_exon_id
where
    seq_region.name in ('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT') and
    ((gene.seq_region_start >= karyotype.seq_region_start and gene.seq_region_start <= karyotype.seq_region_end)
		or (gene.seq_region_end >= karyotype.seq_region_start and gene.seq_region_end <= karyotype.seq_region_end))
	and (t.biotype = 'protein_coding' or display_xref.display_label = 'C11orf95')
    # and display_xref.display_label in ('FAM138A','TP53','CDKN2A','PTEN','MYC')
group by chromosome, gene_start, gene_end, gene_id, gene_name, transcript_id, transcript_version, transcript_start, transcript_end, exon_id, exon_start, exon_end, coding_start, coding_end, strand)
order by if(cast(chromosome as SIGNED) = 0, ascii(chromosome), cast(chromosome as SIGNED)), gene_start, gene_id, exon_start

