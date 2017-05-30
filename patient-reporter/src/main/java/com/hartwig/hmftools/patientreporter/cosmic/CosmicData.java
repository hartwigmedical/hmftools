package com.hartwig.hmftools.patientreporter.cosmic;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CosmicData {
    @NotNull
    private final String entrezId;

    @Nullable
    private final String role;

    @NotNull
    private final String chromosomeBand;

    public CosmicData(@NotNull final String entrezId, @Nullable final String role,
            @NotNull final String chromosomeBand) {
        this.entrezId = entrezId;
        this.role = role;
        this.chromosomeBand = chromosomeBand;
    }

    @NotNull
    public String entrezId() {
        return entrezId;
    }

    @Nullable
    public String role() {
        return role;
    }

    @NotNull
    public String chromosomeBand() {
        return chromosomeBand;
    }
}
