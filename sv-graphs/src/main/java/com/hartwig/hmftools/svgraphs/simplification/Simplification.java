package com.hartwig.hmftools.svgraphs.simplification;

import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;

import com.hartwig.hmftools.svgraphs.BgAdjacency;
import com.hartwig.hmftools.svgraphs.BreakendConsistency;
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
