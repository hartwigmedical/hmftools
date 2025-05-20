package com.hartwig.hmftools.common.collect;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

public class UnionQuickFind<T>
{
    private final HashMap<T, ArrayList<T>> mSetLookup = Maps.newHashMap();
    private final HashMap<T, T> mReprLookup = Maps.newHashMap();

    public void add(final T x)
    {
        if(mSetLookup.containsKey(x))
            return;

        ArrayList<T> set = Lists.newArrayList(List.of(x));
        mSetLookup.put(x, set);
        mReprLookup.put(x, x);
    }

    public T getRepresentative(final T x)
    {
        return mReprLookup.get(x);
    }

    public void merge(final T x, final T y)
    {
        T xRepr = getRepresentative(x);
        T yRepr = getRepresentative(y);
        if(xRepr.equals(yRepr))
            return;

        ArrayList<T> xSet = mSetLookup.get(xRepr);
        ArrayList<T> ySet = mSetLookup.get(yRepr);
        if(xSet.size() >= ySet.size())
        {
            xSet.addAll(ySet);
            mSetLookup.remove(yRepr);
            for(T el : ySet)
                mReprLookup.put(el, xRepr);
        }
        else
        {
            ySet.addAll(xSet);
            mSetLookup.remove(xRepr);
            for(T el : xSet)
                mReprLookup.put(el, yRepr);
        }
    }

    public Collection<T> getSet(final T x)
    {
        return Collections.unmodifiableList(mSetLookup.get(getRepresentative(x)));
    }

    public Collection<Set<T>> getPartitions()
    {
        return mSetLookup.values().stream().map(xs -> (Set<T>) Sets.newHashSet(xs)).toList();
    }
}