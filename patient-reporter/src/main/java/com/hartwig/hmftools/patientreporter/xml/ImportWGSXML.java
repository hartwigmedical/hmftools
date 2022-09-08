package com.hartwig.hmftools.patientreporter.xml;

import java.util.List;

import com.fasterxml.jackson.dataformat.xml.annotation.JacksonXmlElementWrapper;
import com.fasterxml.jackson.dataformat.xml.annotation.JacksonXmlProperty;
import com.hartwig.hmftools.common.xml.KeyXML;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ImportWGSXML {

    @JacksonXmlProperty(localName = "item")
    @JacksonXmlElementWrapper(useWrapping = false)

    @NotNull
    public abstract List<KeyXML> item();
}