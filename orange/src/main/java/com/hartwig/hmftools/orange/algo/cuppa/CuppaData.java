package com.hartwig.hmftools.orange.algo.cuppa;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CuppaData {

    @NotNull
    public abstract List<CuppaPrediction> predictions();

    public abstract int simpleDups32To200B();

    public abstract int maxComplexSize();

    public abstract int telomericSGLs();

    public abstract int LINECount();
}
