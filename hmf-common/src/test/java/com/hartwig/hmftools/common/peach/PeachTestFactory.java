package com.hartwig.hmftools.common.peach;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class PeachTestFactory
{
    @NotNull
    public static ImmutablePeachGenotype.Builder builder()
    {
        return ImmutablePeachGenotype.builder()
                .gene(Strings.EMPTY)
                .allele(Strings.EMPTY)
                .alleleCount(2)
                .function(Strings.EMPTY)
                .linkedDrugs(Strings.EMPTY)
                .urlPrescriptionInfo(Strings.EMPTY)
                .panelVersion(Strings.EMPTY)
                .repoVersion(Strings.EMPTY);
    }
}
