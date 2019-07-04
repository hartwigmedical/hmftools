package com.hartwig.hmftools.linx.visualiser.circos;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositions;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegions;

import org.jetbrains.annotations.NotNull;

public class Span
{

    @NotNull
    public static Map<String, Long> maxPositionPerChromosome(@NotNull final List<GenomePosition> tracks)
    {
        return tracks.stream().collect(Collectors.toMap(GenomePosition::chromosome, GenomePosition::position, Math::max));
    }

    @NotNull
    public static Map<String, Long> minPositionPerChromosome(@NotNull final List<GenomePosition> tracks)
    {
        return tracks.stream().collect(Collectors.toMap(GenomePosition::chromosome, GenomePosition::position, Math::min));
    }

    @NotNull
    public static List<GenomeRegion> span(@NotNull final List<GenomePosition> positions)
    {
        final List<GenomeRegion> result = Lists.newArrayList();

        final List<String> chromosomes = positions.stream().map(GenomePosition::chromosome).distinct().collect(Collectors.toList());
        for (final String chromosome : chromosomes)
        {
            long min =
                    positions.stream().filter(x -> x.chromosome().equals(chromosome)).mapToLong(GenomePosition::position).min().orElse(0);
            long max =
                    positions.stream().filter(x -> x.chromosome().equals(chromosome)).mapToLong(GenomePosition::position).max().orElse(0);

            result.add(GenomeRegions.create(chromosome, min, max));
        }

        Collections.sort(result);
        return result;
    }

    @NotNull
    public static List<GenomeRegion> spanRegions(@NotNull final List<? extends GenomeRegion> regions)
    {
        final List<GenomeRegion> result = Lists.newArrayList();

        final Set<String> chromosomes = regions.stream().map(GenomeRegion::chromosome).collect(Collectors.toSet());
        for (final String chromosome : chromosomes)
        {
            long min = regions.stream().filter(x -> x.chromosome().equals(chromosome)).mapToLong(GenomeRegion::start).min().orElse(0);
            long max = regions.stream().filter(x -> x.chromosome().equals(chromosome)).mapToLong(GenomeRegion::end).max().orElse(0);

            result.add(GenomeRegions.create(chromosome, min, max));
        }

        Collections.sort(result);
        return result;
    }

    @NotNull
    public static List<GenomePosition> allPositions(@NotNull final Collection<? extends GenomeRegion> segments)
    {
        final List<GenomePosition> results = Lists.newArrayList();

        for (final GenomeRegion segment : segments)
        {
            if (HumanChromosome.contains(segment.chromosome()))
            {
                results.add(GenomePositions.create(segment.chromosome(), segment.start()));
                results.add(GenomePositions.create(segment.chromosome(), segment.end()));
            }
        }

        Collections.sort(results);
        return results;
    }

}
