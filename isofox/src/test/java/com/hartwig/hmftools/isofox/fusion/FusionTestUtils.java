package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.isofox.fusion.FusionRead.convertReads;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.isofox.common.ReadRecord;

public final class FusionTestUtils
{
    public static FusionFragment fromReads(final List<ReadRecord> reads)
    {
        return new FusionFragment(new FusionReadGroup(reads.get(0).Id, convertReads(reads)));
    }



}
