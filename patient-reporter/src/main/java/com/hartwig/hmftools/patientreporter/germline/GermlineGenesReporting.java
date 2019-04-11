package com.hartwig.hmftools.patientreporter.germline;

import java.util.Set;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class GermlineGenesReporting {

    @NotNull
    public abstract Set<String> germlineGenes();

    @NotNull
    public abstract Set<String> germlineGenesNotify();

}
