package com.hartwig.hmftools.common.freec;

import com.hartwig.hmftools.common.position.GenomePosition;

import org.immutables.value.Value;

@Value.Style(allParameters = true)
@Value.Immutable
public interface FreecRatio extends GenomePosition {

    double estimatedBAF();

    double ratio();

    double medianRatio();
}
