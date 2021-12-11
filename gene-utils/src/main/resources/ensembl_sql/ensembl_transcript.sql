# Extract transcript and exon data for Ensembl data cache

select q1.*,
if(Strand = -1, ce.seq_region_end - tl.seq_end + 1, cs.seq_region_start + tl.seq_start - 1) as CodingStart,
if(Strand = -1, cs.seq_region_end - tl.seq_start + 1, ce.seq_region_start + tl.seq_end - 1) as CodingEnd
from (
select g.stable_id As GeneId, g.canonical_transcript_id as CanonicalTranscriptId,
t.seq_region_strand as Strand, t.transcript_id as TransId, t.stable_id as Trans, t.biotype as BioType,
t.seq_region_start as TransStart, t.seq_region_end as TransEnd,
et.rank as ExonRank, e.seq_region_start as ExonStart, e.seq_region_end as ExonEnd, e.phase as ExonPhase, e.end_phase as ExonEndPhase
from transcript as t, exon as e, exon_transcript as et, gene as g
where t.transcript_id = et.transcript_id and e.exon_id = et.exon_id
and t.gene_id = g.gene_id
) as q1
left join translation tl on tl.transcript_id = TransId
left join exon cs on cs.exon_id = tl.start_exon_id
left join exon ce on ce.exon_id = tl.end_exon_id
order by GeneId, TransId, ExonStart;
