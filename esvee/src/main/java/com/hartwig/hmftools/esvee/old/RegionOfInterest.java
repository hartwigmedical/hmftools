package com.hartwig.hmftools.esvee.old;

import static java.lang.Math.max;
import static java.lang.Math.min;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.esvee.common.Direction;

public class RegionOfInterest extends ChrBaseRegion
{
    public final Direction Orientation;
    public final boolean Inverted;

    public RegionOfInterest(final String chromosome, final int start, final int end)
    {
        this(chromosome, start, end, Direction.FORWARDS, false);
    }

    public RegionOfInterest(final String chromosome, final int start, final int end, final Direction orientation, final boolean inverted)
    {
        super(chromosome, start, end);
        Orientation = orientation;
        Inverted = inverted;
    }

    public boolean tryExtend(final RegionOfInterest other)
    {
        if(!touches(other))
            return false;

        setStart(min(start(), other.start()));
        setEnd(max(end(), other.end()));
        return true;
    }

    /** Whether this region overlaps or abuts {@code other} */
    public boolean touches(final RegionOfInterest other)
    {
        return Chromosome.equals(other.Chromosome) && Inverted == other.Inverted
                && CommonUtils.touches(start(), end(), other.start(), other.end());
    }

    public int distanceTo(final String chromosome, final int position)
    {
        if(!chromosome.equals(Chromosome))
            return Integer.MAX_VALUE;

        if(start() <= position || position <= end())
            return 0;

        return min(Math.abs(start() - position), Math.abs(end() - position));
    }

    @Override
    public boolean equals(final Object o)
    {
        if(this == o)
            return true;
        if(o == null || getClass() != o.getClass())
            return false;

        final RegionOfInterest region = (RegionOfInterest) o;
        return Chromosome.equals(region.Chromosome)
                && start() == region.start()
                && end() == region.end()
                && Orientation == region.Orientation
                && Inverted == region.Inverted;
    }

    @Override
    public int hashCode()
    {
        int result = Chromosome.hashCode();
        result = 31 * result + start();
        result = 31 * result + end();
        result = 31 * result + Orientation.hashCode();
        result = 31 * result + (Inverted ? 1 : 0);
        return result;
    }

    @Override
    public String toString()
    {
        return String.format("ROI {%s:%d-%d, %s, %s}", Chromosome, start(), end(), Orientation, Inverted);
    }

    public static void appendRegion(final List<RegionOfInterest> existingRegions, final RegionOfInterest region)
    {
        for(RegionOfInterest existingRegion : existingRegions)
            if(existingRegion.tryExtend(region))
                return;

        existingRegions.add(region);
    }

    public static List<RegionOfInterest> tryMerge(final Iterable<RegionOfInterest> regions)
    {
        final Map<String, List<RegionOfInterest>> byChromosome = new LinkedHashMap<>();
        for(RegionOfInterest region : regions)
            byChromosome.computeIfAbsent(region.Chromosome, ignored -> new ArrayList<>()).add(region);

        // Sorting ensures that only a single pass is required to combine the regions, and that we only have to check the last one
        byChromosome.values().forEach(list -> list.sort(Comparator.comparingInt(r -> r.start())));

        return byChromosome.values().stream().flatMap(regionList ->
        {
            final List<RegionOfInterest> consolidated = new ArrayList<>();
            for(RegionOfInterest region : regionList)
                appendRegion(consolidated, region); // TODO: This checks all regions in this chromosome, in theory we only need to check the previous region.
            return consolidated.stream();
        }).collect(Collectors.toList());
    }
}
