package com.hartwig.hmftools.svassembly.models;

public class SimpleSequence implements Sequence
{
    private final String mName;
    private final byte[] mBases;
    private final byte[] mQuals;

    public SimpleSequence(final byte[] bases, final byte[] quals)
    {
        this("Unnamed", bases, quals);
    }

    public SimpleSequence(final String name, final byte[] bases, final byte[] quals)
    {
        mName = name;
        mBases = bases;
        mQuals = quals;
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
    public String toString()
    {
        return getBasesString();
    }
}
