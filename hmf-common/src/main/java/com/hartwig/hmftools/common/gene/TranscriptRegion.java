package com.hartwig.hmftools.common.gene;

import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public interface TranscriptRegion extends GenomeRegion {

    @NotNull
    String transcriptID();

    int transcriptVersion();

    @NotNull
    String gene();

    @NotNull
    String chromosomeBand();

    @NotNull
    default String transcript() {
        return transcriptID() + "." + transcriptVersion();
    }

    @NotNull
    default String chromosomalPosition() { return chromosome() + chromosomeBand(); }
}
