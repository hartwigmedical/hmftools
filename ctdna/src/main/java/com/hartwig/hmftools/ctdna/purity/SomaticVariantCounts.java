package com.hartwig.hmftools.ctdna.purity;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.Integers.median;

import java.util.List;

import com.google.common.collect.Lists;

public class SomaticVariantCounts
{
    public int AlleleFragments;
    public double AllelelQualTotal;

    public final List<Integer> VariantDepths;
    public final List<Integer> NonZeroVariantDepths;

    public SomaticVariantCounts()
    {
        VariantDepths = Lists.newArrayList();
        NonZeroVariantDepths = Lists.newArrayList();
        AlleleFragments = 0;
        AllelelQualTotal = 0;
    }

    public int depthTotal()
    {
        return VariantDepths.stream().mapToInt(x -> x).sum();
    }

    public double medianDepth(boolean useNonZero)
    {
        return useNonZero ? median(NonZeroVariantDepths) : median(VariantDepths);
    }

    public String toString()
    {
        return format("AF(%d) avgQualTotal(%.1f) depthCounts(nonzero=%d all=%d)",
                AlleleFragments, AllelelQualTotal, NonZeroVariantDepths.size(), VariantDepths.size());
    }
}
