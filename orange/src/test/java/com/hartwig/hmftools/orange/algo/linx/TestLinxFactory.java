package com.hartwig.hmftools.orange.algo.linx;

import com.hartwig.hmftools.datamodel.linx.ImmutableLinxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;

public final class TestLinxFactory
{
    public static LinxRecord createMinimalTestLinxData()
    {
        return linxRecordBuilder().build();
    }

    public static ImmutableLinxRecord.Builder linxRecordBuilder()
    {
        return ImmutableLinxRecord.builder();
    }
}
