package com.hartwig.hmftools.datamodel.cuppa;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.util.List;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CuppaData {

    @NotNull
    public abstract List<CuppaPrediction> predictions();

    public abstract int simpleDups32To200B();

    public abstract int maxComplexSize();

    public abstract int telomericSGLs();

    public abstract int lineCount();
}
