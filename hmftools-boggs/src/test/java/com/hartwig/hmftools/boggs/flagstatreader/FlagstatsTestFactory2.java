package com.hartwig.hmftools.boggs.flagstatreader;

import org.jetbrains.annotations.NotNull;

public final class FlagstatsTestFactory2 {

    private FlagstatsTestFactory2() {
    }

    @NotNull
    public static FlagstatData2 createTestData() {
        FlagStatsBuilder builder = new FlagStatsBuilder();
        builder.setMapped(10);
        builder.setTotal(25);
        return new FlagstatData2("AnyPath", builder.build(), builder.build());

    }
}
