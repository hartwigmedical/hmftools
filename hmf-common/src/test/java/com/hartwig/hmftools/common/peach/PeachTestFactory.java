package com.hartwig.hmftools.common.peach;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class PeachTestFactory {

    private PeachTestFactory() {
    }

    @NotNull
    public static ImmutablePeachGenotype.Builder builder() {
        return ImmutablePeachGenotype.builder()
                .gene(Strings.EMPTY)
                .haplotype(Strings.EMPTY)
                .function(Strings.EMPTY)
                .linkedDrugs(Strings.EMPTY)
                .urlPrescriptionInfo(Strings.EMPTY)
                .panelVersion(Strings.EMPTY)
                .repoVersion(Strings.EMPTY);
    }
}
