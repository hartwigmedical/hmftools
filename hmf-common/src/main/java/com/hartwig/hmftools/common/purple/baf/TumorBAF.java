package com.hartwig.hmftools.common.purple.baf;

import com.hartwig.hmftools.common.position.GenomePosition;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface TumorBAF extends GenomePosition {

    double baf();

    default double mBaf() {
        return 0.5 + Math.abs(baf() - 0.5);
    }
}
