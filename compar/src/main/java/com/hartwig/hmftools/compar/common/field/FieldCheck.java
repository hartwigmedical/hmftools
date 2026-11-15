package com.hartwig.hmftools.compar.common.field;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

public class FieldCheck
{
    public final boolean IsCompared;

    public final Double AbsDiffThreshold;
    public final Double PercDiffThreshold;

    // only used for overrides loaded from file
    public final Boolean AbsThresholdOverridden;
    public final Boolean PercThresholdOverridden;

    public FieldCheck(final boolean isCompared)
    {
        this(isCompared, null, null, false, false);
    }

    public FieldCheck(final boolean isCompared, final Double absDiffThreshold, final Double percDiffThreshold)
    {
        this(isCompared, absDiffThreshold, percDiffThreshold, false, false);
    }

    public FieldCheck(
            final boolean isCompared, final Double absDiffThreshold, final Double percDiffThreshold,
            final boolean absThresholdOverridden, final boolean percThresholdOverridden)
    {
        IsCompared = isCompared;
        AbsDiffThreshold = absDiffThreshold;
        PercDiffThreshold = percDiffThreshold;
        AbsThresholdOverridden = absThresholdOverridden;
        PercThresholdOverridden = percThresholdOverridden;
    }

    public boolean hasThresholds() { return AbsDiffThreshold != null || PercDiffThreshold != null; }

    public String toString()
    {
        if(hasThresholds())
            return format("active(%s) thresholds(abs=%.4f, perc=%.4f)", IsCompared, AbsDiffThreshold, PercDiffThreshold);
        else
            return format("active(%s)", IsCompared);
    }

    public boolean hasDifference(double first, double second)
    {
        double absDiff = abs(first - second);
        double percDiff = absDiff / max(abs(first), abs(second));

        boolean satisfiesAbsDiff = AbsDiffThreshold == null || absDiff > AbsDiffThreshold;
        boolean satisfiesRelDiff = PercDiffThreshold == null || percDiff > PercDiffThreshold;

        return satisfiesAbsDiff && satisfiesRelDiff;
    }
}
