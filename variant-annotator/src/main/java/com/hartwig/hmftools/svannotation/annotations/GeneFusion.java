package com.hartwig.hmftools.svannotation.annotations;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class GeneFusion {
    public abstract String type();

    public abstract String start();

    public abstract String geneStart();

    public abstract String geneContextStart();

    public abstract String transcriptStart();

    public abstract String end();

    public abstract String geneEnd();

    public abstract String geneContextEnd();

    public abstract String transcriptEnd();

    public abstract String vaf();
}