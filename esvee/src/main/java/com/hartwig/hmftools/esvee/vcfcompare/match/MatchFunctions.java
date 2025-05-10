package com.hartwig.hmftools.esvee.vcfcompare.match;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;

import com.hartwig.hmftools.esvee.vcfcompare.VariantBreakend;

public class MatchFunctions
{
    public static boolean coordsMatchExactly(final String chrom1, int pos1, byte orient1, final String chrom2, int pos2, byte orient2)
    {
        return orient1 == orient2 && pos1 == pos2 && chrom1.equals(chrom2);
    }

    public static boolean coordsMatchApproximately(
            final String chrom1, int pos1, byte orient1, final String chrom2, int pos2, byte orient2, int upperLowerBounds)
    {
        return orient1 == orient2
            && chrom1.equals(chrom2)
            && positionWithin(pos1, pos2 - upperLowerBounds, pos2 + upperLowerBounds);
    }

    @FunctionalInterface
    public interface MatchFunction
    {
        boolean match(final VariantBreakend breakend1, final VariantBreakend breakend2, boolean checkOtherSide);
    }

    public static class ExactMatcher implements MatchFunction
    {
        @Override
        public boolean match(final VariantBreakend breakend1, final VariantBreakend breakend2, boolean checkOtherSide)
        {
            boolean firstSideMatches =
                    breakend1.Orientation == breakend2.Orientation &&
                    breakend1.minPosition() == breakend2.minPosition() &&
                    breakend1.maxPosition() == breakend2.maxPosition() &&
                    breakend1.Chromosome.equals(breakend2.Chromosome) &&
                    breakend1.InsertSequence.equals(breakend2.InsertSequence) &&
                    breakend1.Homseq.equals(breakend2.Homseq) &&
                    breakend1.SvType.equals(breakend2.SvType);

            if(!checkOtherSide)
                return firstSideMatches;

            boolean secondSideMatches =
                    breakend1.OtherOrientation == breakend2.OtherOrientation &&
                    breakend1.otherMinPosition() == breakend2.otherMinPosition() &&
                    breakend1.otherMaxPosition() == breakend2.otherMaxPosition() &&
                    breakend1.OtherChromosome.equals(breakend2.OtherChromosome);
                    // No need to check insert seq, hom seq or SV type again. It is always the same on the other side

            return firstSideMatches && secondSideMatches;
        }
    }

    public static class CoordsOnlyMatcher implements MatchFunction
    {
        @Override
        public boolean match(final VariantBreakend breakend1, final VariantBreakend breakend2, boolean checkOtherSide)
        {
            boolean firstSideMatches = coordsMatchExactly(
                    breakend1.Chromosome, breakend1.Position, breakend1.Orientation,
                    breakend2.Chromosome, breakend2.Position, breakend2.Orientation
            );

            if(!checkOtherSide)
                return firstSideMatches;

            boolean secondSideMatches = coordsMatchExactly(
                    breakend1.OtherChromosome, breakend1.OtherPosition, breakend1.OtherOrientation,
                    breakend2.OtherChromosome, breakend2.OtherPosition, breakend2.OtherOrientation
            );

            return firstSideMatches && secondSideMatches;
        }
    }

    public static class ApproxMatcher implements MatchFunction
    {
        private static final int UPPER_LOWER_BOUNDS = 10;

        @Override
        public boolean match(final VariantBreakend breakend1, final VariantBreakend breakend2, boolean checkOtherSide)
        {
            boolean firstSideMatches = coordsMatchApproximately(
                    breakend1.Chromosome, breakend1.Position, breakend1.Orientation,
                    breakend2.Chromosome, breakend2.Position, breakend2.Orientation,
                    UPPER_LOWER_BOUNDS
            );

            if(!checkOtherSide)
                return firstSideMatches;

            boolean secondSideMatches = coordsMatchApproximately(
                    breakend1.OtherChromosome, breakend1.OtherPosition, breakend1.OtherOrientation,
                    breakend2.OtherChromosome, breakend2.OtherPosition, breakend2.OtherOrientation,
                    UPPER_LOWER_BOUNDS
            );

            return firstSideMatches && secondSideMatches;

        }
    }
}
