package com.hartwig.hmftools.patientreporter.xml;

import com.fasterxml.jackson.dataformat.xml.annotation.JacksonXmlProperty;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.purple.TumorMutationalStatus;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;

import org.apache.logging.log4j.util.Strings;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class Signature {

    @JacksonXmlProperty(localName = "msscore")
    public abstract double msscore();

    @JacksonXmlProperty(localName = "msstatus")
    @NotNull
    public abstract MicrosatelliteStatus msstatus();

    @JacksonXmlProperty(localName = "tumuload")
    public abstract int tumuload();

    @JacksonXmlProperty(localName = "tumulosta")
    @NotNull
    public abstract TumorMutationalStatus tumulosta();

    @JacksonXmlProperty( localName = "tutmb")
    public abstract double tutmb();

    @JacksonXmlProperty(localName = "horesco")
    public abstract double horesco();

    @JacksonXmlProperty(localName = "horestu")
    @NotNull
    public abstract ChordStatus horestu();

    @JacksonXmlProperty(localName = "geenpv")
    @NotNull
    public abstract String geenpv();


}
