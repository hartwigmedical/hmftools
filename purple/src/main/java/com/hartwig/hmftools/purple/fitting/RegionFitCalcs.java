package com.hartwig.hmftools.purple.fitting;

import static com.hartwig.hmftools.purple.PurpleConstants.MAX_DIPLOID_COPY_NUMBER;
import static com.hartwig.hmftools.purple.PurpleConstants.MIN_DIPLOID_COPY_NUMBER;

import com.hartwig.hmftools.common.utils.Doubles;

public class RegionFitCalcs
{
    public double TumorCopyNumber;
    public double TumorBAF;
    public double RefNormalisedCopyNumber;
    public double MinorAlleleCopyNumberDeviation;
    public double MajorAlleleCopyNumberDeviation;
    public double EventPenalty;
    public double DeviationPenalty;

    public RegionFitCalcs(
            double tumorCopyNumber, double tumorBAF, double refNormalisedCopyNumber, double minorAlleleCopyNumberDeviation,
            double majorAlleleCopyNumberDeviation, double eventPenalty, double deviationPenalty)
    {
        TumorCopyNumber = tumorCopyNumber;
        TumorBAF = tumorBAF;
        RefNormalisedCopyNumber = refNormalisedCopyNumber;
        MinorAlleleCopyNumberDeviation = minorAlleleCopyNumberDeviation;
        MajorAlleleCopyNumberDeviation = majorAlleleCopyNumberDeviation;
        EventPenalty = eventPenalty;
        DeviationPenalty = deviationPenalty;
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
