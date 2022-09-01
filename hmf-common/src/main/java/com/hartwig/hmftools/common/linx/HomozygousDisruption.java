package com.hartwig.hmftools.common.linx;

import com.fasterxml.jackson.dataformat.xml.annotation.JacksonXmlProperty;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class HomozygousDisruption
{
    @JacksonXmlProperty(localName = "chr")
    @NotNull
    public abstract String chromosome();

    @JacksonXmlProperty(localName = "chrbd")
    @NotNull
    public abstract String chromosomeBand();

    @JacksonXmlProperty(localName = "gene")
    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String transcript();

    public abstract boolean isCanonical();
}
