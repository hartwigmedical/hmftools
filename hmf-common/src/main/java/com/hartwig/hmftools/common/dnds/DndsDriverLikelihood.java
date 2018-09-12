package com.hartwig.hmftools.common.dnds;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DndsDriverLikelihood {

    public abstract String gene();

    public abstract double indel();

    public abstract double missense();

    public abstract double nonsense();

    public abstract double splice();

}
