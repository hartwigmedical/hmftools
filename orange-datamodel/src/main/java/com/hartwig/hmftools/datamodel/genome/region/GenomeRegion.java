package com.hartwig.hmftools.datamodel.genome.region;

import org.jetbrains.annotations.NotNull;

public interface GenomeRegion {

    @NotNull
    String chromosome();

    int start();
    int end();
}
