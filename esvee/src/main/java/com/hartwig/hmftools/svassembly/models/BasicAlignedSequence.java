package com.hartwig.hmftools.svassembly.models;

import java.util.List;

public class BasicAlignedSequence implements AlignedSequence
{
    private final String mName;
    private final byte[] mBases;
    private final byte[] mQuals;
    private final List<Alignment> mAlignment;

    public BasicAlignedSequence(final byte[] bases, final byte[] quals, final List<Alignment> alignment)
    {
        this("Unnamed", bases, quals, alignment);
    }

    public BasicAlignedSequence(final String name, final byte[] bases, final byte[] quals, final List<Alignment> alignment)
    {
        mName = name;
        mBases = bases;
        mQuals = quals;
        mAlignment = alignment;
    }

    @Override
    public String getName()
    {
        return mName;
    }

    @Override
    public byte[] getBases()
    {
        return mBases;
    }

    @Override
    public byte[] getBaseQuality()
    {
        return mQuals;
    }

    @Override
    public List<Alignment> getAlignmentBlocks()
    {
        return mAlignment;
    }

    @Override
    public String toString()
    {
        return getBasesString();
    }
}
