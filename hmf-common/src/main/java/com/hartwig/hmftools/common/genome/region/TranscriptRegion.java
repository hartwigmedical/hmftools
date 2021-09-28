package com.hartwig.hmftools.common.genome.region;

import org.jetbrains.annotations.NotNull;

public interface TranscriptRegion extends GenomeRegion {

    @NotNull
    String transName();

    @NotNull
    boolean isCanonical();

    @NotNull
    String geneName();

    @NotNull
    String chromosomeBand();
}
