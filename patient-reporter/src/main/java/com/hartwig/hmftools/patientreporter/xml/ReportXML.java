package com.hartwig.hmftools.patientreporter.xml;

import com.fasterxml.jackson.dataformat.xml.annotation.JacksonXmlProperty;
import com.fasterxml.jackson.dataformat.xml.annotation.JacksonXmlRootElement;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
@JacksonXmlRootElement(localName = "report")
public abstract class ReportXML {

    @JacksonXmlProperty(isAttribute = true, localName = "xmlns:xsi")
    private final String xsi = "http://www.w3.org/2001/XMLSchema-instance";

    @JacksonXmlProperty(isAttribute = true, localName = "xmlns:xsd")
    private final String xsd = "http://www.w3.org/2001/XMLSchema";

    @JacksonXmlProperty(localName = "protocol")
    @NotNull
    public abstract XMLProtocol protocol();
}