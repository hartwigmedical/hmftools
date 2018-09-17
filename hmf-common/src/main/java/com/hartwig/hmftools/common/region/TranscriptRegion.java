package com.hartwig.hmftools.common.region;

import org.jetbrains.annotations.NotNull;

public interface TranscriptRegion extends GenomeRegion {

    @NotNull
    String transcriptID();

    int transcriptVersion();

    @NotNull
    String gene();

    @NotNull
    String chromosomeBand();
}
