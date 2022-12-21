package com.hartwig.hmftools.common.test;

import static java.lang.String.format;

public class ReadIdGenerator
{
    private int mNextReadId;

    public ReadIdGenerator()
    {
        mNextReadId = 1;
    }

    public String nextId()
    {
        String current = currentId();
        mNextReadId++;
        return current;
    }

    public String currentId() { return format("READ_%03d", mNextReadId); }

    public void reset() { mNextReadId = 1;}
}
