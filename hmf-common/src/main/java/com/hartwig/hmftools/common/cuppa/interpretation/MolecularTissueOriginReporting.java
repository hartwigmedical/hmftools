package com.hartwig.hmftools.common.cuppa.interpretation;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularTissueOriginReporting {

    @NotNull
    public abstract String bestCancerType();

    public abstract double bestLikelihood();

    @NotNull
    public abstract String interpretCancerType();

    @Nullable
    public abstract Double interpretLikelihood();

}