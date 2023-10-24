package com.hartwig.hmftools.orange.algo.util;


import com.hartwig.hmftools.datamodel.linx.ImmutableLinxDriver;
import com.hartwig.hmftools.datamodel.linx.LinxDriverType;

import org.jetbrains.annotations.NotNull;

public class LinxDriverTestFactory
{
    @NotNull
    public static ImmutableLinxDriver.Builder builder()
    {
        return ImmutableLinxDriver.builder().gene("").type(LinxDriverType.UNCLEAR);
    }
}
