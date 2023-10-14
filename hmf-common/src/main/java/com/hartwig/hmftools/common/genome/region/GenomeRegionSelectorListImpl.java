package com.hartwig.hmftools.common.genome.region;

import java.util.Collection;
import java.util.List;
import java.util.Optional;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;

import org.jetbrains.annotations.NotNull;

public class GenomeRegionSelectorListImpl<R extends GenomeRegion> implements GenomeRegionSelector<R>
{
    @NotNull
    private final List<R> regions;
    private int index = 0;

    GenomeRegionSelectorListImpl(@NotNull Collection<R> regions)
    {
        this(regions.stream().sorted().collect(Collectors.toList()));
    }

    GenomeRegionSelectorListImpl(@NotNull List<R> regions)
    {
        this.regions = regions;
    }

    @Override
    public void select(@NotNull final GenomeRegion region, @NotNull final Consumer<R> handler)
    {
        if(regions.isEmpty())
        {
            return;
        }

        select(GenomePositions.create(region.chromosome(), region.start()));

        final GenomePosition start = GenomePositions.create(region.chromosome(), region.start());
        final GenomePosition end = GenomePositions.create(region.chromosome(), region.end());

        while(index < regions.size() && compare(start, current()) <= 0 && compare(end, current()) >= 0)
        {
            handler.accept(current());
            index++;
        }

        index = Math.min(index, regions.size() - 1);
    }

    @NotNull
    public Optional<R> select(@NotNull final GenomePosition position)
    {
        if(regions.isEmpty())
        {
            return Optional.empty();
        }

        int currentCompare = compare(position, current());
        while(currentCompare < 0 && index > 0)
        {
            index--;
            currentCompare = compare(position, current());
        }

        while(currentCompare > 0 && index < regions.size() - 1)
        {
            index++;
            currentCompare = compare(position, current());
        }

        return currentCompare == 0 ? Optional.of(current()) : Optional.empty();
    }

    private R current()
    {
        return regions.get(index);
    }

    public static int compare(@NotNull final GenomePosition position, @NotNull final GenomeRegion region)
    {
        int positionChromosome = HumanChromosome.fromString(position.chromosome()).intValue();
        int regionChromosome = HumanChromosome.fromString(region.chromosome()).intValue();
        if(positionChromosome < regionChromosome)
        {
            return -1;
        }
        if(positionChromosome > regionChromosome)
        {
            return 1;
        }

        if(position.position() < region.start())
        {
            return -1;
        }

        if(position.position() > region.end())
        {
            return 1;
        }

        return 0;
    }
}
