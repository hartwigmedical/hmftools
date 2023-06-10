package com.hartwig.hmftools.purple.purity;

import static com.hartwig.hmftools.purple.config.PurpleConstants.MAX_DIPLOID_COPY_NUMBER;
import static com.hartwig.hmftools.purple.config.PurpleConstants.MIN_DIPLOID_COPY_NUMBER;

import com.hartwig.hmftools.common.utils.Doubles;

public class RegionFitCalcs
{
    public final double Purity;
    public final double NormFactor;

    public double TumorCopyNumber;
    public double TumorBAF;
    public double RefNormalisedCopyNumber;
    public double MinorAlleleCopyNumberDeviation;
    public double MajorAlleleCopyNumberDeviation;
    public double EventPenalty;
    public double DeviationPenalty;

    public RegionFitCalcs(final double purity, final double normFactor)
    {
        Purity = purity;
        NormFactor = normFactor;

        TumorCopyNumber = 0;
        TumorBAF = 0;
        RefNormalisedCopyNumber = 0;
        MinorAlleleCopyNumberDeviation = 0;
        MajorAlleleCopyNumberDeviation = 0;
        EventPenalty = 0;
        DeviationPenalty = 0;
    }

    public double majorAlleleCopyNumber() { return TumorBAF * TumorCopyNumber; }
    public double minorAlleleCopyNumber() { return TumorCopyNumber - majorAlleleCopyNumber(); }

    public boolean isDiploid()
    {
        return Doubles.greaterOrEqual(majorAlleleCopyNumber(), MIN_DIPLOID_COPY_NUMBER)
                && Doubles.lessOrEqual(majorAlleleCopyNumber(), MAX_DIPLOID_COPY_NUMBER)
                && Doubles.greaterOrEqual(minorAlleleCopyNumber(), MIN_DIPLOID_COPY_NUMBER)
                && Doubles.lessOrEqual(minorAlleleCopyNumber(), MAX_DIPLOID_COPY_NUMBER);

    }
}
