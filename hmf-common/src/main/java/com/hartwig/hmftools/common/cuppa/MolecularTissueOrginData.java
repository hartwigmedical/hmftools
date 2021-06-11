package com.hartwig.hmftools.common.cuppa;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularTissueOrginData {

    @NotNull
    public abstract String conclusion();

    @NotNull
    public abstract String predictedOrigin();

    @Nullable
    public abstract String predictionLikelihood();
}
