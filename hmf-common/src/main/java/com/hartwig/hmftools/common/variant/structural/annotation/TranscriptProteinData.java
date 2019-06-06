package com.hartwig.hmftools.common.variant.structural.annotation;

public class TranscriptProteinData
{
    public final int TranscriptId;
    public final int TranslationId;
    public final int ProteinFeatureId;
    public final int SeqStart;
    public final int SeqEnd;
    public final String HitDescription;

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
