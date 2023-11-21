package com.hartwig.hmftools.sage.sync;

import htsjdk.samtools.SAMRecord;

public class FragmentData
{
    public final SAMRecord First;
    public final SAMRecord Second;

    public FragmentData(final SAMRecord first, final SAMRecord second)
    {
        First = first;
        Second = second;
    }


}
