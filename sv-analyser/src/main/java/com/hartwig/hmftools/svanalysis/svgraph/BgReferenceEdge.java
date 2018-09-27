package com.hartwig.hmftools.svanalysis.svgraph;

import com.hartwig.hmftools.common.position.GenomePosition;

public class BgReferenceEdge implements BgEdge {
    private final BgSegment left;
    private final BgSegment right;

    public BgReferenceEdge(BgSegment left, BgSegment right) {
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
}
