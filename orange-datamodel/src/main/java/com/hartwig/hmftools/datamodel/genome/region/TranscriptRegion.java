package com.hartwig.hmftools.datamodel.genome.region;

import org.jetbrains.annotations.NotNull;

public interface TranscriptRegion extends GenomeRegion {

    @NotNull
    String transName();

    boolean isCanonical();

    @NotNull
    String geneName();

    @NotNull
    String chromosomeBand();
}
