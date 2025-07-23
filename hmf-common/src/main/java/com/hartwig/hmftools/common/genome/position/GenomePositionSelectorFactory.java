package com.hartwig.hmftools.common.genome.position;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.function.Consumer;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public final class GenomePositionSelectorFactory
{
    private GenomePositionSelectorFactory()
    {
    }

    @NotNull
    public static <P extends GenomePosition> GenomePositionSelector<P> create(@NotNull final Collection<P> positions)
    {
        return new GenomePositionSelectorIteratorImpl<>(positions);
    }

    @NotNull
    public static <P extends GenomePosition> GenomePositionSelector<P> create(@NotNull final List<P> positions)
    {
        return new GenomePositionSelectorListImpl<>(positions);
    }

    @NotNull
    public static <P extends GenomePosition> GenomePositionSelector<P> create(@NotNull final ListMultimap<Chromosome, P> positionsMap)
    {
        return createImpl(new ChromosomePositions<>()
        {
            @Override
            public Set<Chromosome> chromosomes()
            {
                return positionsMap.keySet();
            }

            @Override
            public List<P> positions(final Chromosome chromosome)
            {
                return positionsMap.get(chromosome);
            }
        });
    }

    @NotNull
    public static <P extends GenomePosition> GenomePositionSelector<P> create(@NotNull final Map<Chromosome, List<P>> positions)
    {
        return createImpl(new ChromosomePositions<>()
        {
            @Override
            public Set<Chromosome> chromosomes()
            {
                return positions.keySet();
            }

            @Override
            public List<P> positions(final Chromosome chromosome)
            {
                return positions.get(chromosome);
            }
        });
    }

    @NotNull
    public static <P extends GenomePosition> GenomePositionSelector<P> create(@NotNull final Multimap<Chromosome, P> positions)
    {
        return createImpl(new ChromosomePositions<>()
        {
            @Override
            public Set<Chromosome> chromosomes()
            {
                return positions.keySet();
            }

            @Override
            public GenomePositionSelector<P> createSelector(final Chromosome chromosome)
            {
                return new GenomePositionSelectorIteratorImpl<>(positions.get(chromosome));
            }
        });
    }

    private interface ChromosomePositions<P extends GenomePosition>
    {
        Set<Chromosome> chromosomes();

        default List<P> positions(Chromosome chromosome)
        {
            return Collections.emptyList();
        }

        default GenomePositionSelector<P> createSelector(Chromosome chromosome)
        {
            return new GenomePositionSelectorListImpl<>(positions(chromosome));
        }
    }

    @NotNull
    private static <P extends GenomePosition> GenomePositionSelector<P> createImpl(
            @NotNull final ChromosomePositions<P> chromosomePositions)
    {
        final GenomePositionSelector<P> nullSelector = new NullGenomePositionSelector<>();

        final Map<Chromosome, GenomePositionSelector<P>> chromosomeSelectors = Maps.newHashMap();
        for(final Chromosome chromosome : chromosomePositions.chromosomes())
        {
            chromosomeSelectors.put(chromosome, chromosomePositions.createSelector(chromosome));
        }
        return new GenomePositionSelector<>()
        {
            @NotNull
            @Override
            public Optional<P> select(@NotNull final GenomePosition position)
            {
                return chromosomeSelectors.getOrDefault(position.chr(), nullSelector).select(position);
            }

            @Override
            public void select(final GenomeRegion region, final Consumer<P> handler)
            {
                chromosomeSelectors.getOrDefault(region.chr(), nullSelector).select(region, handler);
            }
        };
    }

    private static class NullGenomePositionSelector<P extends GenomePosition> implements GenomePositionSelector<P>
    {

        @NotNull
        @Override
        public Optional<P> select(@NotNull final GenomePosition position)
        {
            return Optional.empty();
        }

        @Override
        public void select(final GenomeRegion region, final Consumer<P> handler)
        {
            // VOID
        }
    }
}
