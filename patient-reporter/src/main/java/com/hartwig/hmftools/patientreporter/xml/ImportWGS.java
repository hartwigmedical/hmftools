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

    @JacksonXmlProperty(localName = "refNummerWgs")
    @NotNull
    public abstract String refNummerWgs();

    @JacksonXmlProperty(localName = "wgsRedenAanvraag")
    @NotNull
    public abstract String wgsRedenAanvraag();

    @JacksonXmlProperty(localName = "wgsGevrOndzTher")
    @NotNull
    public abstract String wgsGevrOndzTher();

    @JacksonXmlProperty(localName = "wgsGevrOndzTherAnd")
    @NotNull
    public abstract String wgsGevrOndzTherAnd();

    @JacksonXmlProperty(localName = "wgsGevrOndzDiffDiag")
    @NotNull
    public abstract String wgsGevrOndzDiffDiag();

    @JacksonXmlProperty(localName = "wgsGevrOndzDiffDiagAnd")
    @NotNull
    public abstract String wgsGevrOndzDiffDiagAnd();

    @JacksonXmlProperty(localName = "wgsRefNummer")
    @NotNull
    public abstract String wgsRefNummer();

    @JacksonXmlProperty(localName = "wgsPercNeoCellenEx")
    @NotNull
    public abstract String wgsPercNeoCellenEx();

    @JacksonXmlProperty(localName = "wgsPercNeoCellenBeoord")
    @NotNull
    public abstract String wgsPercNeoCellenBeoord();

    @JacksonXmlProperty(localName = "wgsPercNeoCellen")
    @NotNull
    public abstract String wgsPercNeoCellen();

    @JacksonXmlProperty(localName = "wgsDatasheetSeqAnaPanel")
    @NotNull
    public abstract String wgsDatasheetSeqAnaPanel();

    @JacksonXmlProperty(localName = "wgsPlatform")
    @NotNull
    public abstract String wgsPlatform();

    @JacksonXmlProperty(localName = "wgsPlatformAnd")
    @NotNull
    public abstract String wgsPlatformAnd();

    @JacksonXmlProperty(localName = "wgsTumorPurity")
    public abstract double wgsTumorPurity();

    @JacksonXmlProperty(localName = "wgsGemTuPloid")
    public abstract double wgsGemTuPloid();

    @JacksonXmlProperty(localName = "reportableVariants")
    @NotNull
    public abstract List<ReportableVariant> reportableVariants();

    @JacksonXmlProperty(localName = "gainsAndLosses")
    @NotNull
    public abstract List<GainLoss> gainsAndLosses();

    @JacksonXmlProperty(localName = "geneFusions")
    @NotNull
    public abstract List<LinxFusion> geneFusions();

    @JacksonXmlProperty( localName = "signature")
    @NotNull
    public abstract Signature signature();

    @JacksonXmlProperty(localName = "homozygousDisruptions")
    @NotNull
    public abstract List<HomozygousDisruption> homozygousDisruptions();

    @JacksonXmlProperty(localName = "reportableViruses")
    @NotNull
    public abstract List<AnnotatedVirus> reportableViruses();

    @JacksonXmlProperty(localName = "wgsCupAnalyse")
    @NotNull
    public abstract String wgsCupAnalyse();

    @JacksonXmlProperty(localName = "wgsDisclaimerTonen")
    @NotNull
    public abstract String wgsDisclaimerTonen();

    @JacksonXmlProperty(localName = "wgsMolecInter")
    @NotNull
    public abstract String wgsMolecInter();

    @JacksonXmlProperty(localName = "wgsKlinInter")
    @NotNull
    public abstract String wgsKlinInter();

    @JacksonXmlProperty(localName = "wgsAutoKMBP")
    @NotNull
    public abstract String wgsAutoKMBP();

}
