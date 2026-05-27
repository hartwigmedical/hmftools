package com.hartwig.hmftools.common.genome.region;

public class Window
{
    private final int mSize;

    public Window(final int size)
    {
        this.mSize = size;
    }

    public int getSize()
    {
        return mSize;
    }

    // for a standard window of 1000 bases, will return 1001 for values between 1001 - 2000
    public int start(int position)
    {
        return (position - 1) / mSize * mSize + 1;
    }

    // similarly returns 2000 for values 1001 - 2000
    public int end(int position)
    {
        return start(position) + mSize - 1;
    }
}
