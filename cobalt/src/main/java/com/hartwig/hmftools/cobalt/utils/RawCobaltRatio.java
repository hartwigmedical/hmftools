package com.hartwig.hmftools.cobalt.utils;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;

public record RawCobaltRatio(
        @NotNull String chromosome,
        int position,
        double referenceReadCount,
        double tumorReadCount,
        double referenceGcRatio,
        double tumorGcRatio,
        double referenceGcDiploidRatio,
        double referenceGCContent,
        double tumorGCContent
)
{
    public boolean matchesPosition(final RawCobaltRatio other)
    {
        return chromosome.equals(other.chromosome) && position == other.position;
    }

    public RawCobaltRatio differences(final RawCobaltRatio other, final double epsilon)
    {
        Double tumorReadCountDiff = new DoubleDifference(tumorReadCount, other.tumorReadCount, epsilon).difference;
        Double tumorGcRatioDiff = new DoubleDifference(tumorGcRatio, other.tumorGcRatio, epsilon).difference;
        Double tumorGcContentDiff = new DoubleDifference(tumorGCContent, other.tumorGCContent, epsilon).difference;
        if(Doubles.isZero(tumorReadCount))
        {
            Preconditions.checkArgument(Doubles.isZero(other.tumorReadCount));
            tumorReadCountDiff = 0.0;
            tumorGcRatioDiff = 0.0;
            tumorGcContentDiff = 0.0;
        }
        if(Doubles.isZero(tumorGCContent))
        {
            Preconditions.checkArgument(Doubles.isZero(other.tumorGCContent));
            tumorGcRatioDiff = 0.0;
            tumorGcContentDiff = 0.0;
        }
        Double refReadCountDiff = new DoubleDifference(referenceReadCount, other.referenceReadCount, epsilon).difference;
        Double refGcRatioDiff = new DoubleDifference(referenceGcRatio, other.referenceGcRatio, epsilon).difference;
        Double refGcDiploidRatioDiff = new DoubleDifference(referenceGcDiploidRatio, other.referenceGcDiploidRatio, epsilon).difference;
        Double refGcContentDiff = new DoubleDifference(referenceGCContent, other.referenceGCContent, epsilon).difference;
        if(Doubles.isZero(referenceReadCount))
        {
            Preconditions.checkArgument(Doubles.isZero(other.referenceReadCount));
            refReadCountDiff = 0.0;
            refGcRatioDiff = 0.0;
            refGcDiploidRatioDiff = 0.0;
            refGcContentDiff = 0.0;
        }
        if(Doubles.isZero(referenceGCContent))
        {
            Preconditions.checkArgument(Doubles.isZero(other.referenceGCContent));
            refGcRatioDiff = 0.0;
            refGcDiploidRatioDiff = 0.0;
            refGcContentDiff = 0.0;
        }
        boolean hasDifference = !Doubles.isZero(tumorReadCountDiff) ||
                !Doubles.isZero(tumorGcRatioDiff) ||
                !Doubles.isZero(tumorGcContentDiff) ||
                !Doubles.isZero(refReadCountDiff) ||
                !Doubles.isZero(refGcRatioDiff) ||
                !Doubles.isZero(refGcDiploidRatioDiff) ||
                !Doubles.isZero(refGcContentDiff);
        if(hasDifference)
        {
            return new RawCobaltRatio(chromosome, position,
                    refReadCountDiff,
                    tumorReadCountDiff,
                    refGcRatioDiff,
                    tumorGcRatioDiff,
                    refGcDiploidRatioDiff,
                    refGcContentDiff,
                    tumorGcContentDiff);
        }
        return null;
    }
}
