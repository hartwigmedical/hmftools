package com.hartwig.hmftools.boggs.flagstatreader;

import org.jetbrains.annotations.NotNull;

public final class FlagstatsTestFactory {

    private FlagstatsTestFactory() {
    }

    @NotNull
    public static FlagstatData createTestData() {
        FlagStatsBuilder builder = new FlagStatsBuilder();
        builder.setMapped(10);
        builder.setTotal(25);
        return new FlagstatData("AnyPath", builder.build(), builder.build());

    }
}
