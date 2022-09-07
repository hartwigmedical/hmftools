package com.hartwig.hmftools.patientreporter.xml;

import com.fasterxml.jackson.dataformat.xml.annotation.JacksonXmlProperty;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class RubriekXML {

    @JacksonXmlProperty(isAttribute = true)
    private final String naam = "moduleimport";

    @JacksonXmlProperty(localName = "session")
    @NotNull
    public abstract SessionXML jsonSession();
}