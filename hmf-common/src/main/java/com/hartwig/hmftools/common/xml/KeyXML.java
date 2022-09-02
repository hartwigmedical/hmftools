package com.hartwig.hmftools.common.xml;

import java.util.Map;

import com.fasterxml.jackson.dataformat.xml.annotation.JacksonXmlProperty;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class KeyXML {

    @JacksonXmlProperty(localName = "path")
    @NotNull
    public abstract String keyPath();

    @JacksonXmlProperty(localName = "values")
    @NotNull
    public abstract Map<String, String> valuePath();
}