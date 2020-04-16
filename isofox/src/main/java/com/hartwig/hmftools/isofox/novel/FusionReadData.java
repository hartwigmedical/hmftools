package com.hartwig.hmftools.isofox.novel;

import java.util.List;

import com.hartwig.hmftools.isofox.common.ReadRecord;

public class FusionReadData
{
    public final List<ReadRecord> Reads;

    public FusionReadData(List<ReadRecord> reads)
    {
        Reads = reads;
    }

}
