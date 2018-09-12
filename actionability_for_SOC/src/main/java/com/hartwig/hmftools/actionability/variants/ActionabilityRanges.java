package com.hartwig.hmftools.actionability.variants;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })

public abstract class ActionabilityRanges {

    @NotNull
    abstract String gene();

    @NotNull
    abstract String mutationTranscript();

    @NotNull
    abstract String chromosome();

    @NotNull
    abstract String start();

    @NotNull
    abstract String stop();

    @NotNull
    abstract String geneTranscript();

    @NotNull
    abstract String source();

    @NotNull
    abstract String reference();

    @NotNull
    abstract String drugsName();

    @NotNull
    abstract String drugsType();

    @NotNull
    abstract String cancerType();

    @NotNull
    abstract String levelSource();

    @NotNull
    abstract String hmfLevel();

    @NotNull
    abstract String evidenceType();

    @NotNull
    abstract String significanceSource();

    @NotNull
    abstract String hmfResponse();






}
