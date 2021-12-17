package com.hartwig.hmftools.common.genome.position;

import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

public final class GenomePositions {

    private GenomePositions() {
    }

    @NotNull
    public static GenomePosition create(@NotNull final String chromosome, final int position) {
        return ImmutableGenomePositionImpl.builder().chromosome(chromosome).position(position).build();
    }

    @NotNull
    public static <T extends GenomePosition> GenomePosition create(@NotNull final T genomePosition) {
        return ImmutableGenomePositionImpl.builder().chromosome(genomePosition.chromosome()).position(genomePosition.position()).build();
    }

    @NotNull
    public static <S, T extends GenomePosition> Multimap<S, T> union(@NotNull final Multimap<S, T> first,
            @NotNull final Multimap<S, T> second) {
        final Multimap<S, T> union = ArrayListMultimap.create();
        final Set<S> keys = Sets.newHashSet();
        keys.addAll(first.keySet());
        keys.addAll(second.keySet());

        for (S key : keys) {
            final Collection<T> firstCollection = first.get(key);
            final Collection<T> secondCollection = second.get(key);

            if (firstCollection == null) {
                union.putAll(key, secondCollection);
            } else if (secondCollection == null) {
                union.putAll(key, firstCollection);
            } else {
                union.putAll(key, union(firstCollection, secondCollection));
            }
        }

        return union;
    }

    @NotNull
    public static <T extends GenomePosition> List<T> union(@NotNull final Collection<T> first, @NotNull final Collection<T> second) {
        final List<T> merged = Lists.newArrayList();
        final Iterator<T> firstIterator = first.iterator();
        final Iterator<T> secondIterator = second.iterator();

        T firstPosition = firstIterator.hasNext() ? firstIterator.next() : null;
        T secondPosition = secondIterator.hasNext() ? secondIterator.next() : null;
        while (firstPosition != null || secondPosition != null) {

            if (firstPosition == null || (secondPosition != null && secondPosition.compareTo(firstPosition) < 0)) {
                merged.add(secondPosition);
                secondPosition = secondIterator.hasNext() ? secondIterator.next() : null;
            } else if (secondPosition == null || firstPosition.compareTo(secondPosition) < 0) {
                merged.add(firstPosition);
                firstPosition = firstIterator.hasNext() ? firstIterator.next() : null;
            } else {
                merged.add(firstPosition);
                firstPosition = firstIterator.hasNext() ? firstIterator.next() : null;
                secondPosition = secondIterator.hasNext() ? secondIterator.next() : null;
            }
        }

        return merged;
    }
}
