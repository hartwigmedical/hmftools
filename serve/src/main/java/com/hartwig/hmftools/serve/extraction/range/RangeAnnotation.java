package com.hartwig.hmftools.serve.extraction.range;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.serve.extraction.util.MutationTypeFilter;

import org.jetbrains.annotations.NotNull;

public interface RangeAnnotation extends GenomeRegion {

    @NotNull
    String gene();

    @NotNull
    String transcript();

    int rank();

    @NotNull
    MutationTypeFilter mutationType();

}
