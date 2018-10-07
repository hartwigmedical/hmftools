package com.hartwig.hmftools.svanalysis.svgraph;

import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;

public class BgReferenceEdge implements BgEdge {
    private final BgSegment left;
    private final BgSegment right;

    public BgReferenceEdge(BgSegment left, BgSegment right) {
        // TODO: track only positions
        this.left = left;
        this.right = right;
    }

    @Override
    public Double ploidy() {
        return null;
    }

    @Override
    public GenomePosition first() {
        return left;
    }

    @Override
    public GenomePosition second() {
        return right;
    }

    @Override
    public EnrichedStructuralVariant sv() {
        return null;
    }

    @Override
    public String toString() {
        return String.format("REF %s:%d-%d", first().chromosome(), first().position(), second().position());
    }
}
