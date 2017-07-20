package com.hartwig.hmftools.common.purple.ratio;

import com.hartwig.hmftools.common.position.GenomePosition;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface GCContent extends GenomePosition {

    double gcContent();

    double nonNPercentage();

    double mappablePercentage();
}
