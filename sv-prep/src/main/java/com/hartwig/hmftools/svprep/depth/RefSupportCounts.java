package com.hartwig.hmftools.svprep.depth;

class RefSupportCounts
{
    public int RefSupport = 0;
    public int RefPairSupport = 0;

    public RefSupportCounts() {}

    public int total()
    {
        return RefSupport + RefPairSupport;
    }
}
