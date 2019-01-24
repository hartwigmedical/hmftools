# Extract and store all required Ensembl data for VariantAnnotator and SV Analyser
# convert any field with commas to semi-colons, convert tabs to commas


# Extract gene information
# save to ensembl_gene_data.csv

use homo_sapiens_core_89_37;

select gene.stable_id as GeneId,  display_xref.display_label as GeneName, seq_region.name as Chromosome,
	gene.seq_region_strand as Strand, gene.seq_region_start as GeneStart, gene.seq_region_end as GeneEnd,
	GROUP_CONCAT(DISTINCT entrez_xref.dbprimary_acc ORDER BY entrez_xref.dbprimary_acc SEPARATOR ';') as EntrezIds,
    GROUP_CONCAT(DISTINCT karyotype.band ORDER BY karyotype.band SEPARATOR '-') as KaryotypeBand,
	GROUP_CONCAT(DISTINCT syn_xref.dbprimary_acc ORDER BY syn_xref.dbprimary_acc SEPARATOR ';') as Synomyns
from gene
	inner join object_xref as ox on gene.gene_id = ox.ensembl_id and ox.ensembl_object_type = 'GENE'
    inner join xref as display_xref on display_xref.xref_id = gene.display_xref_id
    inner join karyotype on gene.seq_region_id = karyotype.seq_region_id
    inner join seq_region on gene.seq_region_id = seq_region.seq_region_id
    left join xref as entrez_xref on (entrez_xref.xref_id = ox.xref_id and entrez_xref.external_db_id = 1300)
    inner join xref as syn_xref on syn_xref.xref_id = ox.xref_id
where
    seq_region.name in ('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT') and
    ((gene.seq_region_start >= karyotype.seq_region_start and gene.seq_region_start <= karyotype.seq_region_end)
		or (gene.seq_region_end >= karyotype.seq_region_start and gene.seq_region_end <= karyotype.seq_region_end))
	and seq_region.coord_system_id in (select coord_system_id from coord_system where version = 'GRCh37' and name = 'chromosome')
group by Chromosome, GeneStart, GeneEnd, GeneId, GeneName, Strand
order by Chromosome, GeneStart;


# extract transcript and exon information
# save to ensembl_trans_exon_data.csv

select q1.*,
if(Strand = -1, ce.seq_region_end - tl.seq_end + 1, cs.seq_region_start + tl.seq_start - 1) as CodingStart,
if(Strand = -1, cs.seq_region_end - tl.seq_start + 1, ce.seq_region_start + tl.seq_end - 1) as CodingEnd
from (
select g.stable_id As GeneId, g.canonical_transcript_id as CanonicalTranscriptId,
t.seq_region_strand as Strand, t.transcript_id as TransId, t.stable_id as Trans, t.biotype as BioType,
t.seq_region_start as TransStart, t.seq_region_end as TransEnd,
et.rank as ExonRank, e.seq_region_start as ExonStart, e.seq_region_end as ExonEnd, e.phase as ExonPhase, e.end_phase as ExonEndPhase
from transcript as t, exon as e, exon_transcript as et, gene as g, xref as x
where t.transcript_id = et.transcript_id and e.exon_id = et.exon_id and g.display_xref_id = x.xref_id
and t.gene_id = g.gene_id
) as q1
left join translation tl on tl.transcript_id = TransId
left join exon cs on cs.exon_id = tl.start_exon_id
left join exon ce on ce.exon_id = tl.end_exon_id
order by GeneId, TransId, ExonStart;


# extract transcript protein coding information
# save to ensembl_protein_features.csv

select tl.transcript_id, tl.translation_id, protein_feature_id, pf.seq_start, pf.seq_end, hit_description
from protein_feature pf, analysis_description ad, translation tl, transcript t
where pf.analysis_id = ad.analysis_id and pf.translation_id = tl.translation_id and t.transcript_id = tl.transcript_id
and display_label = 'PROSITE profiles'
order by tl.transcript_id, tl.translation_id, pf.seq_start;