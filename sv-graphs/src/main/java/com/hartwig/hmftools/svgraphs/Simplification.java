package com.hartwig.hmftools.svgraphs;

import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;

import org.jetbrains.annotations.NotNull;

import java.util.Collection;
import java.util.List;

public interface Simplification {
    @NotNull
    List<EnrichedStructuralVariant> variants();
    @NotNull
    SimplificationType type();
    @NotNull
    Collection<BreakendConsistency> consistency();
    double copyNumber();
    @NotNull
    List<BgAdjacency> adjacencies();
}
