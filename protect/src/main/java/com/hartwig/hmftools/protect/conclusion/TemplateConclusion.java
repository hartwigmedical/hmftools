package com.hartwig.hmftools.protect.conclusion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class TemplateConclusion {

    @NotNull
    public abstract String abberrationGeneSummary();

    @NotNull
    public abstract String targetedTherapy();

    @NotNull
    public abstract String summaryTextStatement();

    @Nullable
    public abstract String summaryAdditionalText();

}
