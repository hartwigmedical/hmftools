package com.hartwig.hmftools.common.purple.segment;

import java.util.List;

import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Cluster implements GenomeRegion {

    @NotNull
    public abstract List<GenomePosition> ratios();

    @NotNull
    public abstract List<StructuralVariantPosition> variants();
}
