package com.hartwig.hmftools.tars.liftback.features;

import static com.hartwig.hmftools.tars.common.TarsConstants.LOCAL_SV_MAX_LENGTH;

import com.hartwig.hmftools.common.genome.region.Orientation;

// Shared MAPQ-0 preference for discordant mates and supplementary alternatives.
public final class PlacementPairPriority implements Comparable<PlacementPairPriority>
{
    private static final int DELETION = 0;
    private static final int DUPLICATION = 1;
    private static final int INVERSION = 2;
    private static final int OTHER = 3;

    private static final PlacementPairPriority FALLBACK =
            new PlacementPairPriority(OTHER, Long.MAX_VALUE);

    private final int mType;
    private final long mLength;

    private PlacementPairPriority(final int type, final long length)
    {
        mType = type;
        mLength = length;
    }

    static PlacementPairPriority fallback()
    {
        return FALLBACK;
    }

    public static PlacementPairPriority between(
            final String firstChromosome, final int firstPosition, final Orientation firstOrientation,
            final String secondChromosome, final int secondPosition, final Orientation secondOrientation)
    {
        if(!firstChromosome.equals(secondChromosome))
        {
            return FALLBACK;
        }

        long length = Math.abs((long) firstPosition - secondPosition);
        if(length > LOCAL_SV_MAX_LENGTH)
        {
            return FALLBACK;
        }

        if(firstOrientation == secondOrientation)
        {
            return new PlacementPairPriority(INVERSION, length);
        }

        if(length == 0)
        {
            return new PlacementPairPriority(DUPLICATION, length);
        }

        boolean firstIsLower = firstPosition < secondPosition;
        int type = firstIsLower == firstOrientation.isForward() ? DELETION : DUPLICATION;
        return new PlacementPairPriority(type, length);
    }

    @Override
    public int compareTo(final PlacementPairPriority other)
    {
        int typeComparison = Integer.compare(mType, other.mType);
        return typeComparison != 0 ? typeComparison : Long.compare(mLength, other.mLength);
    }
}
