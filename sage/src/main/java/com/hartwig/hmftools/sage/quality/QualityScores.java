package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.sage.quality.QualityCalculator.INVALID_BASE_QUAL;

public class QualityScores
{
    public final double SeqTechBaseQuality;
    public final double RecalibratedBaseQuality;
    public final double FinalBaseQuality;

    public final int FinalMapQuality;

    public final double CombinedQuality;

    public static final QualityScores INVALID_QUAL_SCORES = new QualityScores(
            INVALID_BASE_QUAL, INVALID_BASE_QUAL, 0, INVALID_BASE_QUAL, 0);

    public QualityScores(
            double seqTechBaseQuality, double recalibratedBaseQuality, int finalMapQuality,
            double finalBaseQuality, double combinedQuality)
    {
        SeqTechBaseQuality = seqTechBaseQuality;
        RecalibratedBaseQuality = recalibratedBaseQuality;
        FinalMapQuality = finalMapQuality;
        FinalBaseQuality = finalBaseQuality;
        CombinedQuality = combinedQuality;
    }

    public boolean valid() { return SeqTechBaseQuality != INVALID_BASE_QUAL; }
}
