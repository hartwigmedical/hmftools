package com.hartwig.hmftools.common.collect;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

// see https://en.wikipedia.org/wiki/Disjoint-set_data_structure for an explanation of this data structure
public class UnionFind<T>
{
    private final HashMap<T, T> mParentLookup = Maps.newHashMap();
    private final HashMap<T, Long> mSizeLookup = Maps.newHashMap();

    public void add(final T x)
    {
        if(mParentLookup.containsKey(x))
            return;

        mParentLookup.put(x, null);
        mSizeLookup.put(x, 1L);
    }

    public T getRepresentative(final T x)
    {
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
        T xRep = getRepresentative(x);
        T yRep = getRepresentative(y);
        if(xRep.equals(yRep))
            return;

        long xSize = mSizeLookup.get(xRep);
        long ySize = mSizeLookup.get(yRep);
        long newSize = xSize + ySize;
        if(xSize >= ySize)
        {
            mParentLookup.put(yRep, xRep);
            mSizeLookup.put(xRep, newSize);
        }
        else
        {
            mParentLookup.put(xRep, yRep);
            mSizeLookup.put(yRep, newSize);
        }
    }

    public Collection<Set<T>> getPartitions()
    {
        HashMap<T, Set<T>> partitions = Maps.newHashMap();
        for(T node : mParentLookup.keySet())
        {
            T repr = getRepresentative(node);
            partitions.computeIfAbsent(repr, key -> Sets.newHashSet());
            partitions.get(repr).add(node);
        }

        return partitions.values();
    }
}
