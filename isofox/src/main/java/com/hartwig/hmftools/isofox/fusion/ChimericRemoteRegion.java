package com.hartwig.hmftools.isofox.fusion;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class ChimericRemoteRegion extends ChrBaseRegion
{
    public int Count;

    public ChimericRemoteRegion(final String chromosome, final int regionStart, final int regionEnd)
    {
        super(chromosome, regionStart, regionEnd);
        Count = 1;
    }
}
