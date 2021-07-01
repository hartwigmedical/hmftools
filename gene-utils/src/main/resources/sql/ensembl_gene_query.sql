# Ensemble query is designed to query:
# 1. All protein coding genes.
# 2. All non-protein coding genes that have an entrez id.
# 3. A longer non-canonical transcript of CDKN2A p14arf - a known tumor suppressor.
# 4. C11orf95 - A known fusion candidate in children's ependymoma - added at the request of GARVAN

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
    and gene.stable_id not in
    ('ENSG00000250424','ENSG00000257028','ENSG00000184909','ENSG00000181464','ENSG00000244255','ENSG00000263203','ENSG00000268942','ENSG00000233280',
     'ENSG00000258465','ENSG00000251246','ENSG00000272414','ENSG00000258728','ENSG00000269783','ENSG00000273266','ENSG00000241489','ENSG00000269881',
     'ENSG00000231880','ENSG00000273291','ENSG00000269099','ENSG00000272781','ENSG00000249773','ENSG00000253117','ENSG00000263136','ENSG00000249967',
     'ENSG00000269846','ENSG00000259060','ENSG00000255154','ENSG00000270466','ENSG00000262304','ENSG00000268500','ENSG00000262660','ENSG00000258724',
     'ENSG00000250264','ENSG00000173366','ENSG00000254692','ENSG00000241690','ENSG00000198211','ENSG00000264668','ENSG00000232748','ENSG00000196826',
     'ENSG00000267179','ENSG00000188474','ENSG00000273045','ENSG00000255994','ENSG00000233050','ENSG00000256977','ENSG00000213906','ENSG00000273155',
     'ENSG00000228273','ENSG00000262621','ENSG00000233024','ENSG00000214967','ENSG00000272962','ENSG00000184040','ENSG00000173610','ENSG00000273439')
group by chromosome, gene_start, gene_end, gene_id, gene_name, transcript_id, transcript_version, transcript_start, transcript_end, exon_id, exon_start, exon_end, coding_start, coding_end, strand)
UNION
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
	null as coding_start,
	null as coding_end
from gene
	inner join object_xref on gene.gene_id=object_xref.ensembl_id and object_xref.ensembl_object_type = 'GENE'
    inner join xref as display_xref on (display_xref.xref_id=gene.display_xref_id)
    inner join karyotype on gene.seq_region_id=karyotype.seq_region_id
	inner join transcript t on gene.canonical_transcript_id = t.transcript_id
    inner join seq_region on t.seq_region_id = seq_region.seq_region_id
    inner join exon_transcript et on et.transcript_id = t.transcript_id
    inner join exon e on et.exon_id = e.exon_id
    inner join xref as entrez_xref on (entrez_xref.xref_id=object_xref.xref_id and entrez_xref.external_db_id = 1300)
where
    seq_region.name in ('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT') and
    ((gene.seq_region_start >= karyotype.seq_region_start and gene.seq_region_start <= karyotype.seq_region_end)
		or (gene.seq_region_end >= karyotype.seq_region_start and gene.seq_region_end <= karyotype.seq_region_end))
  and display_xref.display_label not in
  (select display_xref.display_label as gene_name from gene, xref as display_xref, transcript t
	where display_xref.xref_id=gene.display_xref_id and gene.canonical_transcript_id = t.transcript_id and t.biotype = 'protein_coding')
    and t.biotype <> 'protein_coding'
    and display_xref.display_label not in
    ('hsa-mir-3180-3','MIR1587','MIR3615','MIR3916','MIR4461','MIR4519','MIR4523')
group by chromosome, gene_start, gene_end, gene_id, gene_name, transcript_id, transcript_version, transcript_start, transcript_end, exon_id, exon_start, exon_end, coding_start, coding_end, strand)
UNION
(select
 	seq_region.name as chromosome,
     gene.seq_region_start as gene_start,
     gene.seq_region_end as gene_end,
     gene.stable_id as gene_id,
     "CDKN2Ap14ARF" as gene_name,
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
 	inner join transcript t on t.stable_id = 'ENST00000361570'
     inner join seq_region on t.seq_region_id = seq_region.seq_region_id
     inner join exon_transcript et on et.transcript_id = t.transcript_id
     inner join exon e on et.exon_id = e.exon_id
     inner join xref as entrez_xref on (entrez_xref.xref_id=object_xref.xref_id and entrez_xref.external_db_id = 1300)
     inner join translation tl on tl.transcript_id = t.transcript_id
     inner join exon cs on cs.exon_id = tl.start_exon_id
     inner join exon ce on ce.exon_id = tl.end_exon_id
 where
     seq_region.name in ('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT') and
     ((gene.seq_region_start >= karyotype.seq_region_start and gene.seq_region_start <= karyotype.seq_region_end)
 		or (gene.seq_region_end >= karyotype.seq_region_start and gene.seq_region_end <= karyotype.seq_region_end))
 	and t.biotype = 'protein_coding'
     and gene.stable_id = "ENSG00000147889"
group by chromosome, gene_start, gene_end, gene_id, gene_name, transcript_id, transcript_version, transcript_start, transcript_end, exon_id, exon_start, exon_end, coding_start, coding_end, strand)
order by if(cast(chromosome as SIGNED) = 0, ascii(chromosome), cast(chromosome as SIGNED)), gene_start, gene_id, exon_start

