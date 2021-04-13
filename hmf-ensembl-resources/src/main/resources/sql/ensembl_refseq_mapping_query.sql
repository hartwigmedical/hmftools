select
    xref.dbprimary_acc as refSeqId,
    gene.stable_id as geneId,
    transcript.stable_id as transcriptId
from xref
inner join object_xref on xref.xref_id=object_xref.xref_id
inner join transcript on object_xref.ensembl_id=transcript.transcript_id
inner join gene on transcript.gene_id=gene.gene_id
left join xref as entrez_xref on entrez_xref.xref_id=object_xref.xref_id where entrez_xref.external_db_id = 1801;

