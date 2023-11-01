package com.hartwig.hmftools.orange.algo.linx;

import com.hartwig.hmftools.datamodel.linx.ImmutableLinxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;

import org.jetbrains.annotations.NotNull;

public final class TestLinxInterpretationFactory
{
    @NotNull
    public static LinxRecord createMinimalTestLinxData()
    {
        return builder().build();
    }

    @NotNull
    public static ImmutableLinxRecord.Builder builder()
    {
        return ImmutableLinxRecord.builder();
    }
}
