package com.hartwig.hmftools.datamodel.purple;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PurpleFit {

    @NotNull
    public abstract PurpleQC qc();

    public abstract boolean hasSufficientQuality();

    public abstract boolean containsTumorCells();

    public abstract double purity();

    public abstract double ploidy();
}
