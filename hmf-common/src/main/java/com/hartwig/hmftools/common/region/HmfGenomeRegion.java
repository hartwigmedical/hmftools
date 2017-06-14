package com.hartwig.hmftools.common.region;

import org.jetbrains.annotations.NotNull;

public interface HmfGenomeRegion extends GenomeRegion {
    @NotNull
    String transcriptID();

    int transcriptVersion();

    @NotNull
    String gene();

    @NotNull
    String chromosomeBand();

    @NotNull
    String entrezId();

    @NotNull
    default String transcript() {
        return transcriptID() + "." + transcriptVersion();
    }
}
