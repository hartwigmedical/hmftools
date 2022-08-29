package com.hartwig.hmftools.patientreporter.xml;

import java.util.List;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlRootElement;

import com.fasterxml.jackson.dataformat.xml.annotation.JacksonXmlElementWrapper;
import com.fasterxml.jackson.dataformat.xml.annotation.JacksonXmlProperty;
import com.fasterxml.jackson.dataformat.xml.annotation.JacksonXmlRootElement;
import com.hartwig.hmftools.common.linx.HomozygousDisruption;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.purple.loader.GainLoss;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ImportWGS {

    @JacksonXmlProperty(isAttribute = true, localName = "refNummerWgs")
    @NotNull
    public abstract String refNummerWgs();

    @JacksonXmlProperty(isAttribute = true, localName = "wgsRedenAanvraag")
    @NotNull
    public abstract String wgsRedenAanvraag();

    @JacksonXmlProperty(isAttribute = true, localName = "wgsGevrOndzTher")
    @NotNull
    public abstract String wgsGevrOndzTher();

    @JacksonXmlProperty(isAttribute = true, localName = "wgsGevrOndzTherAnd")
    @NotNull
    public abstract String wgsGevrOndzTherAnd();

    @JacksonXmlProperty(isAttribute = true, localName = "wgsGevrOndzDiffDiag")
    @NotNull
    public abstract String wgsGevrOndzDiffDiag();

    @JacksonXmlProperty(isAttribute = true, localName = "wgsGevrOndzDiffDiagAnd")
    @NotNull
    public abstract String wgsGevrOndzDiffDiagAnd();

    @JacksonXmlProperty(isAttribute = true, localName = "wgsRefNummer")
    @NotNull
    public abstract String wgsRefNummer();

    @JacksonXmlProperty(isAttribute = true, localName = "wgsPercNeoCellenEx")
    @NotNull
    public abstract String wgsPercNeoCellenEx();

    @JacksonXmlProperty(isAttribute = true, localName = "wgsPercNeoCellenBeoord")
    @NotNull
    public abstract String wgsPercNeoCellenBeoord();

    @JacksonXmlProperty(isAttribute = true, localName = "wgsPercNeoCellen")
    @NotNull
    public abstract String wgsPercNeoCellen();

    @JacksonXmlProperty(isAttribute = true, localName = "wgsDatasheetSeqAnaPanel")
    @NotNull
    public abstract String wgsDatasheetSeqAnaPanel();

    @JacksonXmlProperty(isAttribute = true, localName = "wgsPlatform")
    @NotNull
    public abstract String wgsPlatform();

    @JacksonXmlProperty(isAttribute = true, localName = "wgsPlatformAnd")
    @NotNull
    public abstract String wgsPlatformAnd();

    @JacksonXmlProperty(isAttribute = true, localName = "wgsTumorPurity")
    public abstract double wgsTumorPurity();

    @JacksonXmlProperty(isAttribute = true, localName = "wgsGemTuPloid")
    public abstract double wgsGemTuPloid();

    @JacksonXmlProperty(isAttribute = true, localName = "reportableVariants")
    @NotNull
    public abstract List<ReportableVariant> reportableVariants();

    @JacksonXmlProperty(isAttribute = true, localName = "gainsAndLosses")
    @NotNull
    public abstract List<GainLoss> gainsAndLosses();

    @JacksonXmlProperty(isAttribute = true, localName = "geneFusions")
    @NotNull
    public abstract List<LinxFusion> geneFusions();

    @JacksonXmlProperty(isAttribute = true, localName = "signature")
    @NotNull
    public abstract Signature signature();

    @JacksonXmlProperty(isAttribute = true, localName = "homozygousDisruptions")
    @NotNull
    public abstract List<HomozygousDisruption> homozygousDisruptions();

    @JacksonXmlProperty(isAttribute = true, localName = "reportableViruses")
    @NotNull
    public abstract List<AnnotatedVirus> reportableViruses();

    @JacksonXmlProperty(isAttribute = true, localName = "wgsCupAnalyse")
    @NotNull
    public abstract String wgsCupAnalyse();

    @JacksonXmlProperty(isAttribute = true, localName = "wgsDisclaimerTonen")
    @NotNull
    public abstract String wgsDisclaimerTonen();

    @JacksonXmlProperty(isAttribute = true, localName = "wgsMolecInter")
    @NotNull
    public abstract String wgsMolecInter();

    @JacksonXmlProperty(isAttribute = true, localName = "wgsKlinInter")
    @NotNull
    public abstract String wgsKlinInter();

    @JacksonXmlProperty(isAttribute = true, localName = "wgsAutoKMBP")
    @NotNull
    public abstract String wgsAutoKMBP();

}
