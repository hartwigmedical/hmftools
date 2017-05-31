package com.hartwig.hmftools.patientreporter.genePanel;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class GenePanelData {
    @NotNull
    private final String entrezId;

    @NotNull
    private final String type;

    @NotNull
    private final String chromosomeBand;

    public GenePanelData(@NotNull final String entrezId, @NotNull final String type,
            @NotNull final String chromosomeBand) {
        this.entrezId = entrezId;
        this.type = type;
        this.chromosomeBand = chromosomeBand;
    }

    @NotNull
    public String entrezId() {
        return entrezId;
    }

    @Nullable
    public String type() {
        return type;
    }

    @NotNull
    public String chromosomeBand() {
        return chromosomeBand;
    }
}
