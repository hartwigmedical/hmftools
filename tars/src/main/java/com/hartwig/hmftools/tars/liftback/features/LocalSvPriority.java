package com.hartwig.hmftools.tars.liftback.features;

import static com.hartwig.hmftools.common.sv.SvUtils.formSvType;
import static com.hartwig.hmftools.tars.common.TarsConstants.LOCAL_SV_MAX_LENGTH;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

// Shared MAPQ-0 placement policy for discordant mates and supplementary alternatives.
public final class LocalSvPriority
{
    private LocalSvPriority() { }

    static Rank fallback()
    {
        return Rank.FALLBACK;
    }

    public static Rank between(
            final String firstChromosome, final int firstPosition, final Orientation firstOrientation,
            final String secondChromosome, final int secondPosition, final Orientation secondOrientation)
    {
        if(!firstChromosome.equals(secondChromosome))
        {
            return Rank.FALLBACK;
        }

        long length = Math.abs((long) firstPosition - secondPosition);
        if(length > LOCAL_SV_MAX_LENGTH)
        {
            return Rank.FALLBACK;
        }

        StructuralVariantType type = formSvType(
                firstChromosome, secondChromosome, firstPosition, secondPosition,
                firstOrientation, secondOrientation, false);
        int typeRank = switch(type)
        {
            case DEL -> 0;
            case DUP -> 1;
            case INV -> 2;
            default -> Rank.FALLBACK_TYPE_RANK;
        };
        return typeRank == Rank.FALLBACK_TYPE_RANK ? Rank.FALLBACK : new Rank(typeRank, length);
    }

    public static final class Rank implements Comparable<Rank>
    {
        private static final int FALLBACK_TYPE_RANK = 3;
        private static final Rank FALLBACK = new Rank(FALLBACK_TYPE_RANK, Long.MAX_VALUE);

        private final int mTypeRank;
        private final long mLength;

        private Rank(final int typeRank, final long length)
        {
            mTypeRank = typeRank;
            mLength = length;
        }

        @Override
        public int compareTo(final Rank other)
        {
            int typeComparison = Integer.compare(mTypeRank, other.mTypeRank);
            return typeComparison != 0 ? typeComparison : Long.compare(mLength, other.mLength);
        }
    }
}
