package com.hartwig.hmftools.orange.algo.linx;

import com.hartwig.hmftools.datamodel.linx.ImmutableLinxDriver;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxDriverEventType;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;

public final class TestLinxFactory
{
    public static LinxRecord createMinimalTestLinxData()
    {
        return linxRecordBuilder().build();
    }

    public static ImmutableLinxDriver.Builder linxDriverBuilder()
    {
        return ImmutableLinxDriver.builder().gene("").type(LinxDriverEventType.UNCLEAR);
    }

    public static ImmutableLinxRecord.Builder linxRecordBuilder()
    {
        return ImmutableLinxRecord.builder();
    }
}
