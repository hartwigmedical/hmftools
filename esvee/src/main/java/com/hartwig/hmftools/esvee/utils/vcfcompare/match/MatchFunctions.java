package com.hartwig.hmftools.esvee.utils.vcfcompare.match;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;

import com.hartwig.hmftools.esvee.utils.vcfcompare.VariantBreakend;

public class MatchFunctions
{
    @FunctionalInterface
    public interface MatchFunction
    {
        boolean match(final VariantBreakend breakend1, final VariantBreakend breakend2);
    }

    public static class ExactMatcher implements MatchFunction
    {
        @Override
        public boolean match(final VariantBreakend breakend1, final VariantBreakend breakend2)
        {
            return
                    // First side
                    breakend1.Orientation == breakend2.Orientation &&
                    breakend1.minPosition() == breakend2.minPosition() &&
                    breakend1.maxPosition() == breakend2.maxPosition() &&
                    breakend1.Chromosome.equals(breakend2.Chromosome) &&
                    breakend1.InsertSequence.equals(breakend2.InsertSequence) &&
                    breakend1.Homseq.equals(breakend2.Homseq) &&
                    breakend1.SvType.equals(breakend2.SvType) &&

                    // Second side
                    breakend1.OtherOrientation == breakend2.OtherOrientation &&
                    breakend1.otherMinPosition() == breakend2.otherMinPosition() &&
                    breakend1.otherMaxPosition() == breakend2.otherMaxPosition() &&
                    breakend1.OtherChromosome.equals(breakend2.OtherChromosome)
                    // No need to check insert seq, hom seq or SV type again. It is always the same on the other side
                    ;
        }
    }

    public static class CoordsOnlyMatcher implements MatchFunction
    {
        @Override
        public boolean match(final VariantBreakend breakend1, final VariantBreakend breakend2)
        {
            return
                    // First side
                    breakend1.Orientation == breakend2.Orientation &&
                    breakend1.Position == breakend2.Position &&
                    breakend1.Chromosome.equals(breakend2.Chromosome) &&


                    // Second side
                    breakend1.OtherOrientation == breakend2.OtherOrientation &&
                    breakend1.OtherPosition == breakend2.OtherPosition &&
                    breakend1.OtherChromosome.equals(breakend2.OtherChromosome)
            ;
        }
    }

    public static class ApproxMatcher implements MatchFunction
    {
        private static final int UPPER_LOWER_BOUNDS = 10;

        @Override
        public boolean match(final VariantBreakend breakend1, final VariantBreakend breakend2)
        {
            return
                    // First side
                    breakend1.Orientation == breakend2.Orientation &&
                    breakend1.Chromosome.equals(breakend2.Chromosome) &&
                    positionWithin(breakend1.Position, breakend2.Position - UPPER_LOWER_BOUNDS, breakend2.Position + UPPER_LOWER_BOUNDS) &&

                    // Second side
                    breakend1.OtherOrientation == breakend2.OtherOrientation &&
                    breakend1.OtherChromosome.equals(breakend2.OtherChromosome) &&
                    positionWithin(breakend1.OtherPosition, breakend2.OtherPosition - UPPER_LOWER_BOUNDS, breakend2.OtherPosition + UPPER_LOWER_BOUNDS)
                    ;
        }
    }
}
