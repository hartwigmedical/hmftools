package com.hartwig.hmftools.bamtools.tofastq;

import static java.lang.String.format;

import htsjdk.samtools.SAMRecord;

public class ReadPair
{
    public final SAMRecord First;
    public final SAMRecord Second;

    public ReadPair(final SAMRecord first, final SAMRecord second)
    {
        First = first;
        Second = second;
    }

    public String toString()
    {
        return format("readId(%s)", First.getReadName());
    }
}
