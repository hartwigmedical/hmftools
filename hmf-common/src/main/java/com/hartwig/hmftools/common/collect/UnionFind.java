package com.hartwig.hmftools.common.collect;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

// see https://en.wikipedia.org/wiki/Disjoint-set_data_structure for an explanation of this data structure
public class UnionFind<T>
{
    private final Map<T, T> mParentLookup;

    public UnionFind()
    {
        mParentLookup = Maps.newHashMap();
    }

    public void add(final T x)
    {
        if(mParentLookup.containsKey(x))
            return;

        mParentLookup.put(x, null);
    }

    public T getRepresentative(final T x)
    {
        add(x);

        T node = x;
        List<T> path = Lists.newArrayList();
        while(node != null)
        {
            path.add(node);
            node = mParentLookup.get(node);
        }

        T repr = path.get(path.size() - 1);
        path.remove(path.size() - 1);
        for(T pathNode : path)
            mParentLookup.put(pathNode, repr);

        return repr;
    }

    public void merge(final T x, final T y)
    {
        add(x);
        add(y);

        T xRep = getRepresentative(x);
        T yRep = getRepresentative(y);
        if(xRep.equals(yRep))
            return;

        mParentLookup.put(yRep, xRep);
    }

    public Collection<Set<T>> getPartitions()
    {
        Map<T, Set<T>> partitions = Maps.newHashMap();
        for(T node : mParentLookup.keySet())
        {
            T repr = getRepresentative(node);
            partitions.computeIfAbsent(repr, key -> Sets.newHashSet());
            partitions.get(repr).add(node);
        }

        return partitions.values();
    }
}
