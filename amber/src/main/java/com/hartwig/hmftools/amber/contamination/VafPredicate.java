package com.hartwig.hmftools.amber.contamination;

interface VafPredicate
{
    boolean test(TumorContamination contamination);
}

class BinomialVafPredicate implements VafPredicate
{
    private final double ContaminationLevel;

    BinomialVafPredicate(double contaminationLevel)
    {
        ContaminationLevel = contaminationLevel;
    }

    @Override
    public boolean test(final TumorContamination contamination)
    {
        final int tumorReadDepth = contamination.Tumor.readDepth();
        final double twoSD = 2 * Math.sqrt(tumorReadDepth * ContaminationLevel * (1 - ContaminationLevel));
        final double mean = ContaminationLevel * tumorReadDepth;
        double rawLowerBound = mean - twoSD;
        double lowerBound = Math.max(Math.ceil(0.003 * tumorReadDepth), rawLowerBound);
        if(contamination.Tumor.altSupport() < lowerBound)
        {
            return false;
        }
        return contamination.Tumor.altSupport() <= mean + twoSD;
    }
}
