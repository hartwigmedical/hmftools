package com.hartwig.hmftools.common.gene;

public class TranscriptProteinData
{
    public final int TranscriptId;
    public final int TranslationId;
    public final int ProteinFeatureId;
    public final int SeqStart;
    public final int SeqEnd;
    public final String HitDescription;

    public static final String BIOTYPE_PROTEIN_CODING = "protein_coding";
    public static final String BIOTYPE_NONSENSE_MED_DECAY = "nonsense_mediated_decay";
    public static final String BIOTYPE_RETAINED_INTRON = "retained_intron";
    public static final String BIOTYPE_PROCESSED_TRANS = "processed_transcript";
    public static final String BIOTYPE_LINC_RNA = "lincRNA";

    public TranscriptProteinData(int transcriptId, int translationId, int proteinFeatureId,
            int seqStart, int seqEnd, final String hitDescription)
    {
        TranscriptId = transcriptId;
        TranslationId = translationId;
        ProteinFeatureId = proteinFeatureId;
        SeqStart = seqStart;
        SeqEnd = seqEnd;
        HitDescription = hitDescription;
    }
}
