package com.hartwig.hmftools.svanalysis.svgraph;

import com.google.common.collect.ImmutableList;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;

import java.util.List;

public class BgSv implements BgEdge {
    private final EnrichedStructuralVariant _sv;

    public BgSv(EnrichedStructuralVariant sv)
    {
        this._sv = sv;
    }

    @Override
    public Double ploidy() {
        return _sv.ploidy();
    }

    @Override
    public GenomePosition first() {
        return _sv.start();
    }

    @Override
    public GenomePosition second() {
        return _sv.end();
    }

    public EnrichedStructuralVariant sv() { return _sv; }

}
