package com.hartwig.hmftools.common.gene;

import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public interface GeneRegion extends GenomeRegion {

    @NotNull
    String transcriptID();

    int transcriptVersion();

    @NotNull
    String gene();

    @NotNull
    default String transcript() {
        return transcriptID() + "." + transcriptVersion();
    }

    @NotNull
    String chromosomeBand();

    @NotNull
    default String chromosomalPosition() { return chromosome() + chromosomeBand(); }
}
