package com.hartwig.hmftools.boggs.flagstatreader;

import org.jetbrains.annotations.NotNull;

public final class FlagStatTestFactory {

    private FlagStatTestFactory() {
    }

    @NotNull
    public static FlagStatData createTestData() {
        FlagStatsBuilder builder = new FlagStatsBuilder();
        builder.setMapped(10);
        builder.setTotal(26);
        return new FlagStatData("AnyPath", builder.build(), builder.build());
    }
}
