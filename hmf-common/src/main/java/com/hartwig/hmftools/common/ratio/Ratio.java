package com.hartwig.hmftools.common.ratio;

import com.hartwig.hmftools.common.position.GenomePosition;
import org.immutables.value.Value;

@Value.Style(allParameters = true)
@Value.Immutable
public interface Ratio extends GenomePosition {

    double ratio();

    double medianRatio();
}
