package com.hartwig.hmftools.svanalysis.svgraph;

import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import org.jetbrains.annotations.NotNull;

import java.util.Collection;
import java.util.List;
import java.util.Set;

public interface Simplification {
    @NotNull
    List<EnrichedStructuralVariant> variants();
    @NotNull
    SimplificationType type();
    @NotNull
    Collection<BreakendConsistency> consistency();
    double ploidy();
    @NotNull
    List<BgAdjacency> adjacencies();
}
