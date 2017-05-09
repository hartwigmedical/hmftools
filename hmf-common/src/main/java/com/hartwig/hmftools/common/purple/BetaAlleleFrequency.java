package com.hartwig.hmftools.common.purple;

import com.hartwig.hmftools.common.position.GenomePosition;
import org.immutables.value.Value;

@Value.Style(allParameters = true)
@Value.Immutable
public abstract class BetaAlleleFrequency implements GenomePosition {

    public abstract double chromosomePosition();

    public abstract double baf();

    public abstract double mBaf();
}
