package com.hartwig.hmftools.amber.contamination;

import com.hartwig.hmftools.amber.VafReading;

interface VafPredicate
{
    boolean test(VafReading contamination);
}

class BinomialVafPredicate implements VafPredicate
{
    private final double ContaminationLevel;

    BinomialVafPredicate(double contaminationLevel)
    {
        ContaminationLevel = contaminationLevel;
    }

    @Override
    public boolean test(final VafReading contamination)
    {
        final int tumorReadDepth = contamination.readDepth();
        final double twoSD = 2 * Math.sqrt(tumorReadDepth * ContaminationLevel * (1 - ContaminationLevel));
        final double mean = ContaminationLevel * tumorReadDepth;
        double rawLowerBound = mean - twoSD;
        double lowerBound = Math.max(Math.ceil(0.003 * tumorReadDepth), rawLowerBound);
        if(contamination.altSupport() < lowerBound)
        {
            return false;
        }
        return contamination.altSupport() <= mean + twoSD;
    }
}
