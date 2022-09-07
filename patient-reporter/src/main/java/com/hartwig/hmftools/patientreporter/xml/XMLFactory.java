package com.hartwig.hmftools.patientreporter.xml;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.linx.HomozygousDisruption;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.purple.loader.CnPerChromosomeArmData;
import com.hartwig.hmftools.common.purple.loader.GainLoss;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.xml.ImmutableKeyItem;
import com.hartwig.hmftools.common.xml.ImmutableKeyXML;
import com.hartwig.hmftools.common.xml.KeyXML;
import com.hartwig.hmftools.patientreporter.algo.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.cfreport.data.GainsAndLosses;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class XMLFactory {

    private XMLFactory() {
    }

    @NotNull
    public static ReportXML generateXMLData(@NotNull AnalysedPatientReport report) {
        ImmutableImportWGSXML importWGS = ImmutableImportWGSXML.builder()
                .refNummerWgs(ImmutableKeyItem.builder()
                        .keyPath(ImmutableKeyXML.builder()
                                .keyPath("refNummerWgs")
                                .valuePath(Map.of("value", report.sampleReport().tumorSampleId()))
                                .build())
                        .build())
                .wgsRedenAanvraag(ImmutableKeyItem.builder()
                        .keyPath(ImmutableKeyXML.builder().keyPath("wgsRedenAanvraag").valuePath(Map.of("value", Strings.EMPTY)).build()).build())
                        .wgsGevrOndzTher(ImmutableKeyXML.builder()
                                .keyPath("wgsGevrOndzTher")
                                .valuePath(Map.of("value", Strings.EMPTY))
                                .build())
                        .wgsGevrOndzTherAnd(ImmutableKeyXML.builder()
                                .keyPath("wgsGevrOndzTherAnd")
                                .valuePath(Map.of("value", Strings.EMPTY))
                                .build())
                        .wgsGevrOndzDiffDiag(ImmutableKeyXML.builder()
                                .keyPath("wgsGevrOndzDiffDiag")
                                .valuePath(Map.of("value", Strings.EMPTY))
                                .build())
                        .wgsGevrOndzDiffDiagAnd(ImmutableKeyXML.builder()
                                .keyPath("wgsGevrOndzDiffDiagAnd")
                                .valuePath(Map.of("value", Strings.EMPTY))
                                .build())
                        .wgsRefNummer(ImmutableKeyXML.builder()
                                .keyPath("wgsRefNummer")
                                .valuePath(Map.of("value", report.sampleReport().tumorSampleId()))
                                .build())
                        .wgsPercNeoCellenEx(ImmutableKeyXML.builder()
                                .keyPath("wgsPercNeoCellenEx")
                                .valuePath(Map.of("value", Strings.EMPTY))
                                .build())
                        .wgsPercNeoCellenBeoord(ImmutableKeyXML.builder()
                                .keyPath("wgsPercNeoCellenBeoord")
                                .valuePath(Map.of("value", Strings.EMPTY))
                                .build())
                        .wgsPercNeoCellen(ImmutableKeyXML.builder()
                                .keyPath("wgsPercNeoCellen")
                                .valuePath(Map.of("value", Strings.EMPTY))
                                .build())
                        .wgsDatasheetSeqAnaPanel(ImmutableKeyXML.builder()
                                .keyPath("wgsDatasheetSeqAnaPanel")
                                .valuePath(Map.of("value", Strings.EMPTY))
                                .build())
                        .wgsPlatform(ImmutableKeyXML.builder().keyPath("wgsPlatform").valuePath(Map.of("value", Strings.EMPTY)).build())
                        .wgsPlatformAnd(ImmutableKeyXML.builder()
                                .keyPath("wgsPlatformAnd")
                                .valuePath(Map.of("value", Strings.EMPTY))
                                .build())
                        .wgsTumorPurity(ImmutableKeyXML.builder()
                                .keyPath("wgsTumorPurity")
                                .valuePath(Map.of("value", Double.toString(report.genomicAnalysis().impliedPurity())))
                                .build())
                        .wgsGemTuPloid(ImmutableKeyXML.builder()
                                .keyPath("wgsGemTuPloid")
                                .valuePath(Map.of("value", Double.toString(report.genomicAnalysis().averageTumorPloidy())))
                                .build())
                        .reportableVariants(addReportableVariantsToXML(report.genomicAnalysis().reportableVariants()))
                        .gainsAndLosses(addGainLossesToXML(report.genomicAnalysis().gainsAndLosses(),
                                report.genomicAnalysis().cnPerChromosome()))
                        .geneFusions(addFusionToXML(report.genomicAnalysis().geneFusions()))
                        .signature(ImmutableSignatureXML.builder()
                                .msscore(ImmutableKeyXML.builder()
                                        .keyPath("importwgs.wgsms.line[1]msscore")
                                        .valuePath(Map.of("value", Double.toString(report.genomicAnalysis().microsatelliteIndelsPerMb())))
                                        .build())
                                .msstatus(ImmutableKeyXML.builder()
                                        .keyPath("importwgs.wgsms.line[1]msstatus")
                                        .valuePath(Map.of("value", report.genomicAnalysis().microsatelliteStatus().toString()))
                                        .build())
                                .tumuload(ImmutableKeyXML.builder()
                                        .keyPath("importwgs.wgsms.line[1]tumuload")
                                        .valuePath(Map.of("value", Double.toString(report.genomicAnalysis().tumorMutationalLoad())))
                                        .build())
                                .tumulosta(ImmutableKeyXML.builder()
                                        .keyPath("importwgs.wgsms.line[1]tumulosta")
                                        .valuePath(Map.of("value", report.genomicAnalysis().tumorMutationalLoadStatus().toString()))
                                        .build())
                                .tutmb(ImmutableKeyXML.builder()
                                        .keyPath("importwgs.wgsms.line[1]tutmb")
                                        .valuePath(Map.of("value", Double.toString(report.genomicAnalysis().tumorMutationalBurden())))
                                        .build())
                                .horesco(ImmutableKeyXML.builder()
                                        .keyPath("importwgs.wgsms.line[1]horesco")
                                        .valuePath(Map.of("value", Double.toString(report.genomicAnalysis().chordHrdValue())))
                                        .build())
                                .horestu(ImmutableKeyXML.builder()
                                        .keyPath("importwgs.wgsms.line[1]horestu")
                                        .valuePath(Map.of("value", report.genomicAnalysis().chordHrdStatus().toString()))
                                        .build())
                                .geenpv(ImmutableKeyXML.builder()
                                        .keyPath("importwgs.wgsms.line[1]geenpv")
                                        .valuePath(Map.of("value", Strings.EMPTY))
                                        .build())
                                .build())
                        .homozygousDisruptions(addHomozygousDisruptionsToXML(report.genomicAnalysis().homozygousDisruptions()))
                        .reportableViruses(addVirussesToXML(report.genomicAnalysis().reportableViruses()))
                        .wgsCupAnalyse(ImmutableKeyXML.builder().keyPath("wgsCupAnalyse").valuePath(Map.of("value", Strings.EMPTY)).build())
                        .wgsDisclaimerTonen(ImmutableKeyXML.builder()
                                .keyPath("wgsDisclaimerTonen")
                                .valuePath(Map.of("value", Strings.EMPTY))
                                .build())
                        .wgsMolecInter(ImmutableKeyXML.builder().keyPath("wgsMolecInter").valuePath(Map.of("value", Strings.EMPTY)).build())
                        .wgsKlinInter(ImmutableKeyXML.builder().keyPath("wgsKlinInter").valuePath(Map.of("value", Strings.EMPTY)).build())
                        .wgsAutoKMBP(ImmutableKeyXML.builder().keyPath("wgsAutoKMBP").valuePath(Map.of("value", Strings.EMPTY)).build())
                        .build();

        return ImmutableReportXML.builder()
                .protocol(ImmutableXMLProtocol.builder()
                        .meta(ImmutableProtocolNameXML.builder().protocolName("Moleculairebepalingen").build())
                        .content(ImmutableContentXML.builder()
                                .rubriek(ImmutableRubriekXML.builder()
                                        .jsonSession(ImmutableSessionXML.builder().importWGSNew(importWGS).build())
                                        .build())
                                .build())
                        .build())
                .build();
    }

    @NotNull
    public static List<KeyXML> addGainLossesToXML(@NotNull List<GainLoss> gainLosses,
            @NotNull List<CnPerChromosomeArmData> chromosomeArmData) {
        List<KeyXML> keyXMLList = Lists.newArrayList();
        int count = 1;
        for (GainLoss gainLoss : gainLosses) {
            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgscnv.line[" + count + "]export")
                    .valuePath(Map.of("value", Strings.EMPTY))
                    .build());

            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgscnv.line[" + count + "]chr")
                    .valuePath(Map.of("value", gainLoss.chromosome()))
                    .build());

            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgscnv.line[" + count + "]region")
                    .valuePath(Map.of("value", gainLoss.chromosomeBand()))
                    .build());

            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgscnv.line[" + count + "]gene")
                    .valuePath(Map.of("value", gainLoss.gene()))
                    .build());
            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgscnv.line[" + count + "]type")
                    .valuePath(Map.of("value", gainLoss.interpretation().display()))
                    .build());

            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgscnv.line[" + count + "]copies")
                    .valuePath(Map.of("value", Double.toString(gainLoss.minCopies())))
                    .build());

            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgscnv.line[" + count + "]charmco")
                    .valuePath(Map.of("value", GainsAndLosses.chromosomeArmCopyNumber(chromosomeArmData, gainLoss)))
                    .build());

            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgscnv.line[" + count + "]geenpv")
                    .valuePath(Map.of("value", Strings.EMPTY))
                    .build());

            count += 1;
        }
        return keyXMLList;
    }

    @NotNull
    public static List<KeyXML> addHomozygousDisruptionsToXML(@NotNull List<HomozygousDisruption> homozygousDisruptions) {
        List<KeyXML> keyXMLList = Lists.newArrayList();
        int count = 1;
        for (HomozygousDisruption homozygousDisruption : homozygousDisruptions) {
            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importw.wgsgshzy.line[" + count + "]export")
                    .valuePath(Map.of("value", Strings.EMPTY))
                    .build());

            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgshzy.line[" + count + "]gen")
                    .valuePath(Map.of("value", homozygousDisruption.gene()))
                    .build());

            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgshzy.line[" + count + "]chr")
                    .valuePath(Map.of("value", homozygousDisruption.chromosome()))
                    .build());

            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgshzy.line[" + count + "]chrbd")
                    .valuePath(Map.of("value", homozygousDisruption.chromosomeBand()))
                    .build());

            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgshzy.line[" + count + "]geenpv")
                    .valuePath(Map.of("value", Strings.EMPTY))
                    .build());

            count += 1;
        }
        return keyXMLList;
    }

    @NotNull
    public static List<KeyXML> addReportableVariantsToXML(@NotNull List<ReportableVariant> reportableVariants) {
        List<KeyXML> keyXMLList = Lists.newArrayList();
        int count = 1;
        for (ReportableVariant reportableVariant : reportableVariants) {
            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgsgene.line[" + count + "]export")
                    .valuePath(Map.of("value", Strings.EMPTY))
                    .build());

            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgsgene.line[" + count + "]name")
                    .valuePath(Map.of("value", reportableVariant.gene()))
                    .build());

            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgsgene.line[" + count + "]pos")
                    .valuePath(Map.of("value", reportableVariant.gDNA()))
                    .build());

            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgsgene.line[" + count + "]var")
                    .valuePath(Map.of("value", reportableVariant.canonicalHgvsCodingImpact()))
                    .build());

            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgsgene.line[" + count + "]prot")
                    .valuePath(Map.of("value", reportableVariant.canonicalHgvsProteinImpact()))
                    .build());

            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgsgene.line[" + count + "]readep")
                    .valuePath(Map.of("value", reportableVariant.alleleReadCount() + "/" + reportableVariant.totalReadCount()))
                    .build());

            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgsgene.line[" + count + "]copie")
                    .valuePath(Map.of("value", Double.toString(reportableVariant.totalCopyNumber())))
                    .build());

            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgsgene.line[" + count + "]tvaf")
                    .valuePath(Map.of("value", reportableVariant.tVAF()))
                    .build());

            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgsgene.line[" + count + "]biallc")
                    .valuePath(Map.of("value",
                            reportableVariant.biallelic() == null ? Strings.EMPTY : Boolean.toString(reportableVariant.biallelic())))
                    .build());

            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgsgene.line[" + count + "]hotsp")
                    .valuePath(Map.of("value", reportableVariant.hotspot().name()))
                    .build());

            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgsgene.line[" + count + "]driver")
                    .valuePath(Map.of("value", reportableVariant.driverLikelihoodInterpretation().display()))
                    .build());

            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgsgene.line[" + count + "]geenpv")
                    .valuePath(Map.of("value", Strings.EMPTY))
                    .build());
            count += 1;
        }
        return keyXMLList;
    }

    @NotNull
    public static List<KeyXML> addVirussesToXML(@NotNull List<AnnotatedVirus> annotatedVirusList) {
        List<KeyXML> keyXMLList = Lists.newArrayList();
        int count = 1;
        for (AnnotatedVirus virus : annotatedVirusList) {
            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgsvrs.line[" + count + "]export")
                    .valuePath(Map.of("value", Strings.EMPTY))
                    .build());

            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgsvrs.line[" + count + "]name")
                    .valuePath(Map.of("value", virus.name()))
                    .build());

            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgsvrs.line[" + count + "]geenpv")
                    .valuePath(Map.of("value", Strings.EMPTY))
                    .build());
            count += 1;
        }
        return keyXMLList;
    }

    @NotNull
    public static List<KeyXML> addFusionToXML(@NotNull List<LinxFusion> linxFusions) {
        List<KeyXML> keyXMLList = Lists.newArrayList();
        int count = 1;
        for (LinxFusion fusion : linxFusions) {
            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgsfusie.line[" + count + "]export")
                    .valuePath(Map.of("value", Strings.EMPTY))
                    .build());

            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgsfusie.line[" + count + "]name")
                    .valuePath(Map.of("value", fusion.name()))
                    .build());
            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgsfusie.line[" + count + "]f5gen")
                    .valuePath(Map.of("value", fusion.geneStart()))
                    .build());
            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgsfusie.line[" + count + "]f5refid")
                    .valuePath(Map.of("value", fusion.geneTranscriptStart()))
                    .build());
            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgsfusie.line[" + count + "]f5exon")
                    .valuePath(Map.of("value", fusion.geneContextStart()))
                    .build());
            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgsfusie.line[" + count + "]f3gen")
                    .valuePath(Map.of("value", fusion.geneEnd()))
                    .build());
            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgsfusie.line[" + count + "]f3refid")
                    .valuePath(Map.of("value", fusion.geneTranscriptEnd()))
                    .build());
            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgsfusie.line[" + count + "]f3exon")
                    .valuePath(Map.of("value", fusion.geneContextEnd()))
                    .build());
            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgsfusie.line[" + count + "]tuco")
                    .valuePath(Map.of("value", Double.toString(fusion.junctionCopyNumber())))
                    .build());
            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgsfusie.line[" + count + "]fufra")
                    .valuePath(Map.of("value", fusion.phased().displayStr()))
                    .build());
            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgsfusie.line[" + count + "]driver")
                    .valuePath(Map.of("value", fusion.likelihood().displayStr()))
                    .build());
            keyXMLList.add(ImmutableKeyXML.builder()
                    .keyPath("importwgs.wgsfusie.line[" + count + "]geenpv")
                    .valuePath(Map.of("value", Strings.EMPTY))
                    .build());
            count += 1;
        }
        return keyXMLList;
    }
}