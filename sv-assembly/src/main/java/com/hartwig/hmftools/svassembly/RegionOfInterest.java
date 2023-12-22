package com.hartwig.hmftools.svassembly;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import com.hartwig.hmftools.common.sv.Direction;
import com.hartwig.hmftools.svassembly.util.RangeUtils;

public class RegionOfInterest
{
    public final String Chromosome;
    public int Start;
    public int End;

    public final Direction Orientation;
    public final boolean Inverted;

    public RegionOfInterest(final String chromosome, final int start, final int end)
    {
        this(chromosome, start, end, Direction.FORWARDS, false);
    }

    public RegionOfInterest(final String chromosome, final int start, final int end, final Direction orientation, final boolean inverted)
    {
        assert start <= end : "Start must not be after end for RegionOfInterest";
        Chromosome = chromosome;
        Orientation = orientation;
        Inverted = inverted;
        Start = start;
        End = end;
    }

    public boolean tryExtend(final RegionOfInterest other)
    {
        if(!touches(other))
            return false;

        Start = Math.min(Start, other.Start);
        End = Math.max(End, other.End);
        return true;
    }

    /** Whether this region overlaps or abuts {@code other} */
    public boolean touches(final RegionOfInterest other)
    {
        return Chromosome.equals(other.Chromosome) && Inverted == other.Inverted
                && RangeUtils.touches(Start, End, other.Start, other.End);
    }

    public int distanceTo(final String chromosome, final int position)
    {
        if (!chromosome.equals(Chromosome))
            return Integer.MAX_VALUE;

        if (Start <= position || position <= End)
            return 0;

        return Math.min(Math.abs(Start - position), Math.abs(End - position));
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
                && Start == region.Start
                && End == region.End
                && Orientation == region.Orientation
                && Inverted == region.Inverted;
    }

    @Override
    public int hashCode()
    {
        int result = Chromosome.hashCode();
        result = 31 * result + Start;
        result = 31 * result + End;
        result = 31 * result + Orientation.hashCode();
        result = 31 * result + (Inverted ? 1 : 0);
        return result;
    }

    @Override
    public String toString()
    {
        return String.format("ROI {%s:%d-%d, %s, %s}", Chromosome, Start, End, Orientation, Inverted);
    }

    public static void appendRegion(final List<RegionOfInterest> existingRegions, final RegionOfInterest region)
    {
        for(final RegionOfInterest existingRegion : existingRegions)
            if(existingRegion.tryExtend(region))
                return;

        existingRegions.add(region);
    }

    public static List<RegionOfInterest> tryMerge(final List<? extends Iterable<RegionOfInterest>> regions)
    {
        return tryMerge(regions.stream().flatMap(r -> StreamSupport.stream(r.spliterator(), false))::iterator);
    }

    public static List<RegionOfInterest> tryMerge(final Iterable<RegionOfInterest> regions)
    {
        final Map<String, List<RegionOfInterest>> byChromosome = new LinkedHashMap<>();
        for(final RegionOfInterest region : regions)
            byChromosome.computeIfAbsent(region.Chromosome, ignored -> new ArrayList<>()).add(region);

        // Sorting ensures that only a single pass is required to combine the regions, and that we only have to check the last one
        byChromosome.values().forEach(list -> list.sort(Comparator.comparingInt(r -> r.Start)));

        return byChromosome.values().stream().flatMap(regionList ->
        {
            final List<RegionOfInterest> consolidated = new ArrayList<>();
            for(final RegionOfInterest region : regionList)
                appendRegion(consolidated, region); // TODO: This checks all regions in this chromosome, in theory we only need to check the previous region.
            return consolidated.stream();
        }).collect(Collectors.toList());
    }
}
