package com.hartwig.hmftools.actionability.compare_with_SOC;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })

public abstract class ReadingMethodsFile {

    @NotNull
    abstract String cancerType();

    @NotNull
    abstract String LMS();

    @NotNull
    abstract String detectingDeviation();

    @NotNull
    abstract String technic();

    @NotNull
    abstract String moleculairTest();

    @NotNull
    abstract String percentageTumoCel();

    @NotNull
    abstract String results();
}
