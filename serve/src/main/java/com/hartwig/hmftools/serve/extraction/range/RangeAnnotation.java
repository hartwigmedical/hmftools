package com.hartwig.hmftools.serve.extraction.range;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.serve.extraction.util.MutationTypeFilter;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public interface RangeAnnotation extends GenomeRegion {

    @NotNull
    String gene();

    @NotNull
    String transcript();

    @Nullable
    Integer rank();

    @NotNull
    MutationTypeFilter mutationType();

}
