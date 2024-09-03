package com.hartwig.hmftools.common.metrics;

import org.jetbrains.annotations.NotNull;

public final class FlagstatTestFactory
{
    @NotNull
    public static BamFlagStats createMinimalTestFlagstat()
    {
        return ImmutableBamFlagStats.builder()
                .uniqueReadCount(0)
                .secondaryCount(0)
                .supplementaryCount(0)
                .duplicateProportion(0D)
                .mappedProportion(0D)
                .pairedInSequencingProportion(0D)
                .properlyPairedProportion(0D)
                .withItselfAndMateMappedProportion(0D)
                .singletonProportion(0D)
                .build();
    }
}
