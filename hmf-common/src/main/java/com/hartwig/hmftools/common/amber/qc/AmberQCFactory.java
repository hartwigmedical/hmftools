package com.hartwig.hmftools.common.amber.qc;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class AmberQCFactory {

    private AmberQCFactory() {
    }

    @NotNull
    public static AmberQC create(double contamination, double consanguinityProportion, @Nullable String uniparentalDisomy) {

        return ImmutableAmberQC.builder()
                .contamination(contamination)
                .consanguinityProportion(consanguinityProportion)
                .uniparentalDisomy(uniparentalDisomy).build();
    }
}
