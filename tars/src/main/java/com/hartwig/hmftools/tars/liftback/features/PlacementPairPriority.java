package com.hartwig.hmftools.tars.liftback.features;

import static com.hartwig.hmftools.tars.common.TarsConstants.LOCAL_SV_MAX_LENGTH;

import com.hartwig.hmftools.common.genome.region.Orientation;

// Shared MAPQ-0 preference for discordant mates and supplementary alternatives.
public final class PlacementPairPriority implements Comparable<PlacementPairPriority>
{
    // mates this close are one fragment rather than a rearrangement, whatever their orientations imply
    private static final long LOCAL_PAIR_MAX_SEPARATION = 1000;

    private static final int LOCAL = 0;
    private static final int DELETION = 1;
    private static final int DUPLICATION = 2;
    private static final int INVERSION = 3;
    private static final int OTHER = 4;

    private static final PlacementPairPriority FALLBACK =
            new PlacementPairPriority(OTHER, Long.MAX_VALUE);

    private final int mCategory;
    private final long mLength;

    private PlacementPairPriority(final int category, final long length)
    {
        mCategory = category;
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
        return between(
                firstChromosome, firstPosition, firstOrientation,
                secondChromosome, secondPosition, secondOrientation,
                Math.abs((long) firstPosition - secondPosition));
    }

    // only for mate pairs: a supplementary alternative sits this close to its primary by construction
    public static PlacementPairPriority betweenMates(
            final String firstChromosome, final int firstPosition, final Orientation firstOrientation,
            final String secondChromosome, final int secondPosition, final Orientation secondOrientation,
            final long separation)
    {
        if(firstChromosome.equals(secondChromosome) && separation <= LOCAL_PAIR_MAX_SEPARATION)
        {
            return new PlacementPairPriority(LOCAL, separation);
        }

        return between(
                firstChromosome, firstPosition, firstOrientation, secondChromosome, secondPosition, secondOrientation,
                separation);
    }

    // separation ranks the pair; the breakend positions still decide its category
    public static PlacementPairPriority between(
            final String firstChromosome, final int firstPosition, final Orientation firstOrientation,
            final String secondChromosome, final int secondPosition, final Orientation secondOrientation,
            final long separation)
    {
        if(!firstChromosome.equals(secondChromosome))
        {
            return FALLBACK;
        }

        if(separation > LOCAL_SV_MAX_LENGTH)
        {
            return FALLBACK;
        }

        if(firstOrientation == secondOrientation)
        {
            return new PlacementPairPriority(INVERSION, separation);
        }

        long breakendLength = Math.abs((long) firstPosition - secondPosition);
        if(breakendLength == 0)
        {
            return new PlacementPairPriority(DUPLICATION, separation);
        }

        boolean firstIsLower = firstPosition < secondPosition;
        int category = firstIsLower == firstOrientation.isForward() ? DELETION : DUPLICATION;
        return new PlacementPairPriority(category, separation);
    }

    @Override
    public int compareTo(final PlacementPairPriority other)
    {
        int categoryComparison = Integer.compare(mCategory, other.mCategory);
        return categoryComparison != 0 ? categoryComparison : Long.compare(mLength, other.mLength);
    }
}
