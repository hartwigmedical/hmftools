package com.hartwig.hmftools.redux.common;

import htsjdk.samtools.SAMRecord;

public class ReadWithFragCoords
{
    public final SAMRecord Read;
    public final FragmentCoords FragCoords;

    public ReadWithFragCoords(final SAMRecord read, final FragmentCoords fragCoords)
    {
        Read = read;
        FragCoords = fragCoords;
    }
}
