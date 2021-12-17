package com.hartwig.hmftools.common.utils.pcf;

import com.hartwig.hmftools.common.genome.position.GenomePosition;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PCFPosition implements GenomePosition {

    @NotNull
    public abstract PCFSource source();

    public abstract int minPosition();

    public abstract int maxPosition();
}
