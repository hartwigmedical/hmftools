package com.hartwig.hmftools.purple.fitting;

import com.hartwig.hmftools.purple.region.ObservedRegion;

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
        return ObservedRegion.isDiploid(majorAlleleCopyNumber()) && ObservedRegion.isDiploid(minorAlleleCopyNumber());

    }
}
