package com.hartwig.hmftools.common.purple;

import com.hartwig.hmftools.common.region.GenomeRegion;

import org.immutables.value.Value;

@Value.Immutable
public abstract class ConsolidatedRegion implements GenomeRegion {

    public abstract double averageBAF();

    public abstract double averageRatioOfRatios();
}
