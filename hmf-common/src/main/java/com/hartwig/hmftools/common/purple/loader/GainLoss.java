package com.hartwig.hmftools.common.purple.loader;

import com.fasterxml.jackson.dataformat.xml.annotation.JacksonXmlProperty;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class GainLoss {

    @JacksonXmlProperty(isAttribute = true, localName = "type")
    @NotNull
    public abstract CopyNumberInterpretation interpretation();

    @JacksonXmlProperty(isAttribute = true, localName = "chr")
    @NotNull
    public abstract String chromosome();

    @JacksonXmlProperty(isAttribute = true, localName = "region")
    @NotNull
    public abstract String chromosomeBand();

    @JacksonXmlProperty(isAttribute = true, localName = "gene")
    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String transcript();

    public abstract boolean isCanonical();

    @JacksonXmlProperty(isAttribute = true, localName = "minCopies")
    public abstract long minCopies();

    @JacksonXmlProperty(isAttribute = true, localName = "maxCopies")
    public abstract long maxCopies();
}
