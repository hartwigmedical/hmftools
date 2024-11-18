package com.hartwig.hmftools.fastqtools.biomodalcollapse;

import java.util.StringJoiner;

public class StringJoinerCounter
{
    private final StringJoiner mStringJoiner;

    private int mComponentCount;

    public StringJoinerCounter(final String delimiter)
    {
        mStringJoiner = new StringJoiner(delimiter);
        mComponentCount = 0;
    }

    public void add(final String s)
    {
        mStringJoiner.add(s);
        mComponentCount++;
    }

    public int componentCount()
    {
        return mComponentCount;
    }

    @Override
    public String toString()
    {
        return mStringJoiner.toString();
    }
}
