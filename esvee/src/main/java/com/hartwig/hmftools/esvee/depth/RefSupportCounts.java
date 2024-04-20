package com.hartwig.hmftools.esvee.depth;

import static java.lang.String.format;

public class RefSupportCounts
{
    public int RefSupport = 0;
    public int RefPairSupport = 0;
    public final int VafCap;

    public RefSupportCounts(final int vafCap) { VafCap = vafCap; }

    public int total()
    {
        return RefSupport + RefPairSupport;
    }

    public boolean exceedsMaxDepth() { return VafCap > 0 && RefSupport + RefPairSupport >= VafCap; }

    public String toString()
    {
        return format("support(ref=%d pair=%d) vafCap(%d)", RefSupport, RefPairSupport, VafCap);
    }
}
