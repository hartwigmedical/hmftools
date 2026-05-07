package com.hartwig.hmftools.common.genome.position;

import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

public final class GenomePositions
{
    public static GenomePosition create(final String chromosome, final int position)
    {
        return new GenomePositionImpl(chromosome, position);
    }

    public static <T extends GenomePosition> GenomePosition create(final T genomePosition)
    {
        return new GenomePositionImpl(genomePosition.chromosome(), genomePosition.position());
    }

    public static <S, T extends GenomePosition> Multimap<S, T> union(final Multimap<S, T> first, final Multimap<S, T> second)
    {
        final Multimap<S, T> union = ArrayListMultimap.create();
        final Set<S> keys = Sets.newHashSet();
        keys.addAll(first.keySet());
        keys.addAll(second.keySet());

        for(S key : keys)
        {
            final Collection<T> firstCollection = first.get(key);
            final Collection<T> secondCollection = second.get(key);

            if(firstCollection == null)
            {
                union.putAll(key, secondCollection);
            }
            else if(secondCollection == null)
            {
                union.putAll(key, firstCollection);
            }
            else
            {
                union.putAll(key, union(firstCollection, secondCollection));
            }
        }

        return union;
    }

    public static <S, T extends GenomePosition> Map<S,List<T>> unionOfMaps(final Map<S,List<T>> first, final Map<S,List<T>> second)
    {
        Map<S,List<T>> unionMap = Maps.newHashMap();
        final Set<S> keys = Sets.newHashSet();
        keys.addAll(first.keySet());
        keys.addAll(second.keySet());

        for(S key : keys)
        {
            List<T> firstCollection = first.get(key);
            List<T> secondCollection = second.get(key);

            if(firstCollection == null)
            {
                unionMap.put(key, secondCollection);
            }
            else if(secondCollection == null)
            {
                unionMap.put(key, firstCollection);
            }
            else
            {
                unionMap.put(key, union(firstCollection, secondCollection));
            }
        }

        return unionMap;
    }

    public static <T extends GenomePosition> List<T> union(final Collection<T> first, final Collection<T> second)
    {
        final List<T> merged = Lists.newArrayList();
        final Iterator<T> firstIterator = first.iterator();
        final Iterator<T> secondIterator = second.iterator();

        T firstPosition = firstIterator.hasNext() ? firstIterator.next() : null;
        T secondPosition = secondIterator.hasNext() ? secondIterator.next() : null;
        while(firstPosition != null || secondPosition != null)
        {

            if(firstPosition == null || (secondPosition != null && secondPosition.compareTo(firstPosition) < 0))
            {
                merged.add(secondPosition);
                secondPosition = secondIterator.hasNext() ? secondIterator.next() : null;
            }
            else if(secondPosition == null || firstPosition.compareTo(secondPosition) < 0)
            {
                merged.add(firstPosition);
                firstPosition = firstIterator.hasNext() ? firstIterator.next() : null;
            }
            else
            {
                merged.add(firstPosition);
                firstPosition = firstIterator.hasNext() ? firstIterator.next() : null;
                secondPosition = secondIterator.hasNext() ? secondIterator.next() : null;
            }
        }

        return merged;
    }
}
