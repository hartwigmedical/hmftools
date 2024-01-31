package com.hartwig.hmftools.esvee.common;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.esvee.read.Read;

public class SequenceMismatches
{
    private Map<Integer,BaseMismatches> mBaseMismatches;

    public SequenceMismatches()
    {
        mBaseMismatches = null;
    }

    public boolean hasMismatches() { return mBaseMismatches != null; }

    public int positionCount() { return mBaseMismatches != null ? mBaseMismatches.size() : 0; }

    public int distinctBaseCount()
    {
        return mBaseMismatches != null ? mBaseMismatches.values().stream().mapToInt(x -> x.mismatchCount()).sum() : 0;
    }

    public int readBaseCount()
    {
        return mBaseMismatches != null ? mBaseMismatches.values().stream().mapToInt(x -> x.mismatchReadTotal()).sum() : 0;
    }

    public void add(final int assemblyIndex, final byte base, final Read read, final byte qual)
    {
        if(mBaseMismatches == null)
        {
            mBaseMismatches = Maps.newLinkedHashMap();
            BaseMismatches baseMismatches = new BaseMismatches(new BaseMismatch(base, read, qual));
            mBaseMismatches.put(assemblyIndex, baseMismatches);
            return;
        }

        BaseMismatches existing = mBaseMismatches.get(assemblyIndex);

        if(existing != null)
        {
            existing.addMismatch(base, read, qual);
        }
        else
        {
            BaseMismatches baseMismatches = new BaseMismatches(new BaseMismatch(base, read, qual));
            mBaseMismatches.put(assemblyIndex, baseMismatches);
        }
    }

    public Map<Integer,BaseMismatches> indexedBaseMismatches() { return mBaseMismatches; }

    @VisibleForTesting
    public List<BaseMismatch> allBaseMismatches()
    {
        if(mBaseMismatches == null)
            return Collections.emptyList();

        List<BaseMismatch> baseMismatches = Lists.newArrayList();

        mBaseMismatches.values().forEach(x -> baseMismatches.addAll(x.baseMismatches()));
        return baseMismatches;
    }
}
