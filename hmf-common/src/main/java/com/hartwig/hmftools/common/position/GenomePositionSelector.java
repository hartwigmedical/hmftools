package com.hartwig.hmftools.common.position;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.region.GenomeRegion;

public interface GenomePositionSelector<P extends GenomePosition>  {
    void select(GenomeRegion region, Consumer<P> handler);
}
