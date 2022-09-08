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

    //    @JacksonXmlProperty(isAttribute = true)
    //    @NotNull
    //    public abstract KeyItem refNummerWgs();
    //
    //    @JacksonXmlProperty(isAttribute = true)
    //    @NotNull
    //    public abstract KeyItem wgsRedenAanvraag();
    //
    //    @JacksonXmlProperty(localName = "itemWgsGevrOndzTher")
    //    @NotNull
    //    public abstract KeyXML wgsGevrOndzTher();
    //
    //    @JacksonXmlProperty(localName = "itemWgsGevrOndzTherAnd")
    //    @NotNull
    //    public abstract KeyXML wgsGevrOndzTherAnd();
    //
    //    @JacksonXmlProperty(localName = "itemWgsGevrOndzDiffDiag")
    //    @NotNull
    //    public abstract KeyXML wgsGevrOndzDiffDiag();
    //
    //    @JacksonXmlProperty(localName = "itemWgsGevrOndzDiffDiagAnd")
    //    @NotNull
    //    public abstract KeyXML wgsGevrOndzDiffDiagAnd();
    //
    //    @JacksonXmlProperty(localName = "itemWgsRefNummer")
    //    @NotNull
    //    public abstract KeyXML wgsRefNummer();
    //
    //    @JacksonXmlProperty(localName = "itemWgsPercNeoCellenEx")
    //    @NotNull
    //    public abstract KeyXML wgsPercNeoCellenEx();
    //
    //    @JacksonXmlProperty(localName = "itemWgsPercNeoCellenBeoord")
    //    @NotNull
    //    public abstract KeyXML wgsPercNeoCellenBeoord();
    //
    //    @JacksonXmlProperty(localName = "itemWgsPercNeoCellen")
    //    @NotNull
    //    public abstract KeyXML wgsPercNeoCellen();
    //
    //    @JacksonXmlProperty(localName = "itemWgsDatasheetSeqAnaPanel")
    //    @NotNull
    //    public abstract KeyXML wgsDatasheetSeqAnaPanel();
    //
    //    @JacksonXmlProperty(localName = "itemWgsPlatform")
    //    @NotNull
    //    public abstract KeyXML wgsPlatform();
    //
    //    @JacksonXmlProperty(localName = "itemWgsPlatformAnd")
    //    @NotNull
    //    public abstract KeyXML wgsPlatformAnd();
    //
    //    @JacksonXmlProperty(localName = "itemWgsTumorPurity")
    //    public abstract KeyXML wgsTumorPurity();
    //
    //    @JacksonXmlProperty(localName = "itemWgsGemTuPloid")
    //    public abstract KeyXML wgsGemTuPloid();
    //
    //    @JacksonXmlProperty(localName = "itemReportableVariants")
    //    @JacksonXmlElementWrapper(useWrapping = false)
    //    @NotNull
    //    public abstract List<KeyXML> reportableVariants();
    //
    //    @JacksonXmlProperty(localName = "itemGainsAndLosses")
    //    @JacksonXmlElementWrapper(useWrapping = false)
    //    @NotNull
    //    public abstract List<KeyXML> gainsAndLosses();
    //
    //    @JacksonXmlProperty(localName = "itemGeneFusions")
    //    @JacksonXmlElementWrapper(useWrapping = false)
    //    @NotNull
    //    public abstract List<KeyXML> geneFusions();
    //
    //    @JacksonXmlProperty(localName = "itemSignature")
    //    @NotNull
    //    public abstract SignatureXML signature();
    //
    //    @JacksonXmlProperty(localName = "itemHomozygousDisruptions")
    //    @JacksonXmlElementWrapper(useWrapping = false)
    //    @NotNull
    //    public abstract List<KeyXML> homozygousDisruptions();
    //
    //    @JacksonXmlProperty(localName = "itemReportableViruses")
    //    @JacksonXmlElementWrapper(useWrapping = false)
    //    @NotNull
    //    public abstract List<KeyXML> reportableViruses();
    //
    //    @JacksonXmlProperty(localName = "itemWgsCupAnalyse")
    //    @NotNull
    //    public abstract KeyXML wgsCupAnalyse();
    //
    //    @JacksonXmlProperty(localName = "itemWgsDisclaimerTonen")
    //    @NotNull
    //    public abstract KeyXML wgsDisclaimerTonen();
    //
    //    @JacksonXmlProperty(localName = "itemWgsMolecInter")
    //    @NotNull
    //    public abstract KeyXML wgsMolecInter();
    //
    //    @JacksonXmlProperty(localName = "itemWgsKlinInter")
    //    @NotNull
    //    public abstract KeyXML wgsKlinInter();
    //
    //    @JacksonXmlProperty(localName = "itemWgsAutoKMBP")
    //    @NotNull
    //    public abstract KeyXML wgsAutoKMBP();
}