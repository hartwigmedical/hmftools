package com.hartwig.hmftools.sage.seqtech;

import java.util.List;
import java.util.Objects;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

@VisibleForTesting
public class MergedHomopolymers
{
    public final List<Homopolymer> RefHomopolymers;
    public final List<Homopolymer> ReadHomopolymers;

    private boolean mVariantInMergedHomopolymers;
    private List<RefMask> mRefMasks;

    public MergedHomopolymers(final List<Homopolymer> refHomopolymers, final List<Homopolymer> readHomopolymers,
            boolean variantInMergedHomopolymers)
    {
        RefHomopolymers = refHomopolymers;
        ReadHomopolymers = readHomopolymers;
        mVariantInMergedHomopolymers = variantInMergedHomopolymers;
        mRefMasks = Lists.newArrayList();
    }

    public MergedHomopolymers()
    {
        this(Lists.newArrayList(), Lists.newArrayList(), false);
    }

    public void setVariantInMergedHomopolymers()
    {
        mVariantInMergedHomopolymers = true;
    }

    public boolean variantInMergedHomopolymers()
    {
        return mVariantInMergedHomopolymers;
    }

    public void addRefMask(final RefMask refMask)
    {
        mRefMasks.add(refMask);
    }

    public List<RefMask> refMasks()
    {
        return mRefMasks;
    }

    // TODO: remove
    @Override
    public boolean equals(final Object o)
    {
        if(this == o)
            return true;

        if(!(o instanceof MergedHomopolymers))
            return false;

        final MergedHomopolymers that = (MergedHomopolymers) o;
        return mVariantInMergedHomopolymers == that.mVariantInMergedHomopolymers
                && Objects.equals(RefHomopolymers, that.RefHomopolymers)
                && Objects.equals(ReadHomopolymers, that.ReadHomopolymers) && Objects.equals(mRefMasks, that.mRefMasks);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(RefHomopolymers, ReadHomopolymers, mVariantInMergedHomopolymers, mRefMasks);
    }
}
