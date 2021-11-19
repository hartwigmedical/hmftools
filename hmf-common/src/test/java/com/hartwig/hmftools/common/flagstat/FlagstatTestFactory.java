package com.hartwig.hmftools.common.flagstat;

import org.jetbrains.annotations.NotNull;

public final class FlagstatTestFactory {

    private FlagstatTestFactory() {
    }

    @NotNull
    public static Flagstat createMinimalTestFlagstat() {
        return ImmutableFlagstat.builder()
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
