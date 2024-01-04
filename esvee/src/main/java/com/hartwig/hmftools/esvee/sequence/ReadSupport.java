package com.hartwig.hmftools.esvee.sequence;

import com.hartwig.hmftools.esvee.read.Read;

public class ReadSupport
{
    public final Read Read;
    public final int Index;

    public ReadSupport(final com.hartwig.hmftools.esvee.read.Read read, final int index)
    {
        Read = read;
        Index = index;
    }
}
