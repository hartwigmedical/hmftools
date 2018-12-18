package com.hartwig.hmftools.common.actionability;

import com.hartwig.hmftools.common.actionability.ActionabilitySource;
import com.hartwig.hmftools.common.actionability.EvidenceScope;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ClinicalTrial {

    @NotNull
    public abstract String event();

    @NotNull
    public abstract String acronym();

    @NotNull
    public abstract ActionabilitySource source();

    @NotNull
    public abstract String reference();

    public abstract boolean isOnLabel();

    @NotNull
    public abstract String cancerType();

    @NotNull
    public abstract EvidenceScope scope();
    
    @NotNull
    public abstract String type();

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String chromosome();

    @NotNull
    public abstract String position();

    @NotNull
    public abstract String ref();

    @NotNull
    public abstract String alt();

    @NotNull
    public abstract String cnvType();

    @NotNull
    public abstract String fusionFiveGene();

    @NotNull
    public abstract String fusionThreeGene();

}
