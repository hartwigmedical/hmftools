# Extract transcript protein data for Ensembl data cache

select tl.transcript_id as TranscriptId, tl.translation_id as TranslationId, protein_feature_id as ProteinFeatureId,
pf.seq_start as SeqStart, pf.seq_end as SeqEnd, hit_description as HitDescription
from protein_feature pf, analysis_description ad, translation tl, transcript t
where pf.analysis_id = ad.analysis_id and pf.translation_id = tl.translation_id and t.transcript_id = tl.transcript_id
and display_label = 'PROSITE profiles'
order by tl.transcript_id, tl.translation_id, pf.seq_start;
