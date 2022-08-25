package com.hartwig.hmftools.common.peach;

import com.fasterxml.jackson.dataformat.xml.annotation.JacksonXmlProperty;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PeachGenotype {

    @JacksonXmlProperty(isAttribute = true, localName = "gene")
    @NotNull
    public abstract String gene();

    @JacksonXmlProperty(isAttribute = true, localName = "haplotype")
    @NotNull
    public abstract String haplotype();

    @JacksonXmlProperty(isAttribute = true, localName = "function")
    @NotNull
    public abstract String function();

    @JacksonXmlProperty(isAttribute = true, localName = "linkedDrugs")
    @NotNull
    public abstract String linkedDrugs();

    @JacksonXmlProperty(isAttribute = true, localName = "urlPrescriptionInfo")
    @NotNull
    public abstract String urlPrescriptionInfo();

    @JacksonXmlProperty(isAttribute = true, localName = "panelVersion")
    @NotNull
    public abstract String panelVersion();

    @JacksonXmlProperty(isAttribute = true, localName = "repoVersion")
    @NotNull
    public abstract String repoVersion();

}
