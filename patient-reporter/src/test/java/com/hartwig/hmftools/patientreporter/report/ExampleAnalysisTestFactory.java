package com.hartwig.hmftools.patientreporter.report;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestFactory.createTestCopyNumberBuilder;
import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testBaseReportData;
import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testSampleReport;

import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.actionability.ActionabilitySource;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.EvidenceLevel;
import com.hartwig.hmftools.common.actionability.ImmutableEvidenceItem;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.ImmutableFittedPurity;
import com.hartwig.hmftools.common.variant.Clonality;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.BaseReportData;
import com.hartwig.hmftools.patientreporter.ImmutableAnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.actionability.ClinicalTrial;
import com.hartwig.hmftools.patientreporter.actionability.ImmutableClinicalTrial;
import com.hartwig.hmftools.patientreporter.chord.ChordAnalysis;
import com.hartwig.hmftools.patientreporter.chord.ImmutableChordAnalysis;
import com.hartwig.hmftools.patientreporter.germline.GermlineVariant;
import com.hartwig.hmftools.patientreporter.structural.ImmutableReportableGeneDisruption;
import com.hartwig.hmftools.patientreporter.structural.ReportableGeneDisruption;
import com.hartwig.hmftools.patientreporter.structural.ReportableGeneFusion;
import com.hartwig.hmftools.patientreporter.variants.ImmutableReportableSomaticVariant;
import com.hartwig.hmftools.patientreporter.variants.ReportableSomaticVariant;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

final class ExampleAnalysisTestFactory {

    private ExampleAnalysisTestFactory() {
    }

    @NotNull
    static AnalysedPatientReport buildCOLO829() {
        final double pathologyTumorPercentage = 0.8;
        final double impliedTumorPurity = 1D;
        final double averageTumorPloidy = 3.1;
        final int tumorMutationalLoad = 182;
        final double tumorMutationalBurden = 13.6;
        final double microsatelliteIndelsPerMb = 0.1089;

        final BaseReportData baseReportData = testBaseReportData();
        final FittedPurity fittedPurity = createFittedPurity(impliedTumorPurity, averageTumorPloidy);

        final PurityAdjuster purityAdjuster = new PurityAdjuster(Gender.MALE, fittedPurity);
        final List<EvidenceItem> clinicalEvidence = createCOLO829ClinicalEvidence();
        final List<ClinicalTrial> clinicalTrials = createCOLO829ClinicalTrials();
        final List<ReportableSomaticVariant> somaticVariants = createCOLO829SomaticVariants(purityAdjuster);
        final List<GeneCopyNumber> copyNumbers = createCOLO829CopyNumbers();
        final List<ReportableGeneFusion> fusions = createCOLO829Fusions();
        final List<GermlineVariant> germlineVariants = createCOLO829GermlineVariants(purityAdjuster);
        final ChordAnalysis chordAnalysis = createCOLO829ChordAnalysis();
        final List<ReportableGeneDisruption> disruptions = createCOLO829Disruptions();

        final SampleReport sampleReport = testSampleReport(pathologyTumorPercentage);

        return ImmutableAnalysedPatientReport.of(sampleReport,
                fittedPurity.purity(),
                true,
                fittedPurity.ploidy(),
                clinicalEvidence,
                clinicalTrials,
                somaticVariants,
                microsatelliteIndelsPerMb,
                tumorMutationalLoad,
                tumorMutationalBurden,
                chordAnalysis,
                true,
                germlineVariants,
                copyNumbers,
                fusions,
                disruptions,
                Resources.getResource("circos/circos_example.png").getPath(),
                Optional.of("this is a test report and does not relate to any real patient"),
                baseReportData.signaturePath(),
                baseReportData.logoRVAPath());
    }

    @NotNull
    private static List<EvidenceItem> createCOLO829ClinicalEvidence() {
        List<EvidenceItem> evidenceItems = Lists.newArrayList();
        evidenceItems.add(evidenceBuilder().event("BRAF p.Val600Glu")
                .drug("Dabrafenib")
                .level(EvidenceLevel.LEVEL_A)
                .response("Responsive")
                .reference("V600E")
                .source(ActionabilitySource.ONCOKB)
                .build());

        evidenceItems.add(evidenceBuilder().event("BRAF p.Val600Glu")
                .drug("Dabrafenib + Trametinib")
                .level(EvidenceLevel.LEVEL_A)
                .response("Responsive")
                .reference("BRAF:V600E,V600K")
                .source(ActionabilitySource.CGI)
                .build());

        evidenceItems.add(evidenceBuilder().event("BRAF p.Val600Glu")
                .drug("Encorafenib + Binimetinib")
                .level(EvidenceLevel.LEVEL_A)
                .response("Responsive")
                .reference("V600E")
                .source(ActionabilitySource.ONCOKB)
                .build());

        evidenceItems.add(evidenceBuilder().event("BRAF p.Val600Glu")
                .drug("Trametinib")
                .level(EvidenceLevel.LEVEL_A)
                .response("Responsive")
                .reference("BRAF:V600E,V600K")
                .source(ActionabilitySource.CGI)
                .build());

        evidenceItems.add(evidenceBuilder().event("BRAF p.Val600Glu")
                .drug("Vemurafenib")
                .level(EvidenceLevel.LEVEL_A)
                .response("Responsive")
                .reference("BRAF:V600E,V600D,V600K,V600M,V600G,V600R")
                .source(ActionabilitySource.CGI)
                .build());

        evidenceItems.add(evidenceBuilder().event("BRAF p.Val600Glu")
                .drug("Vemurafenib + Cobimetinib")
                .level(EvidenceLevel.LEVEL_A)
                .response("Responsive")
                .reference("variant:17")
                .source(ActionabilitySource.CIVIC)
                .build());

        evidenceItems.add(evidenceBuilder().event("BRAF p.Val600Glu")
                .drug("Alpelisib + Cetuximab + Encorafenib")
                .level(EvidenceLevel.LEVEL_B)
                .response("Responsive")
                .reference("variant:17")
                .source(ActionabilitySource.CIVIC)
                .build());

        evidenceItems.add(evidenceBuilder().event("BRAF p.Val600Glu")
                .drug("Bevacizumab")
                .level(EvidenceLevel.LEVEL_B)
                .response("Resistant")
                .reference("variant:12")
                .source(ActionabilitySource.CIVIC)
                .build());

        evidenceItems.add(evidenceBuilder().event("BRAF p.Val600Glu")
                .drug("Cetuximab")
                .level(EvidenceLevel.LEVEL_B)
                .response("Resistant")
                .reference("BRAF:V600E")
                .source(ActionabilitySource.CGI)
                .build());

        evidenceItems.add(evidenceBuilder().event("BRAF p.Val600Glu")
                .drug("Encorafenib + Cetuximab")
                .level(EvidenceLevel.LEVEL_B)
                .response("Responsive")
                .reference("variant:17")
                .source(ActionabilitySource.CIVIC)
                .build());

        evidenceItems.add(evidenceBuilder().event("BRAF p.Val600Glu")
                .drug("Irinotecan")
                .level(EvidenceLevel.LEVEL_B)
                .response("Resistant")
                .reference("variant:399")
                .source(ActionabilitySource.CIVIC)
                .build());

        evidenceItems.add(evidenceBuilder().event("BRAF p.Val600Glu")
                .drug("Oxaliplatin")
                .level(EvidenceLevel.LEVEL_B)
                .response("Resistant")
                .reference("variant:12")
                .source(ActionabilitySource.CIVIC)
                .build());

        evidenceItems.add(evidenceBuilder().event("BRAF p.Val600Glu")
                .drug("Panitumumab")
                .level(EvidenceLevel.LEVEL_B)
                .response("Resistant")
                .reference("variant:17")
                .source(ActionabilitySource.CIVIC)
                .build());

        evidenceItems.add(evidenceBuilder().event("BRAF p.Val600Glu")
                .drug("RO4987655")
                .level(EvidenceLevel.LEVEL_B)
                .response("Responsive")
                .reference("variant:17")
                .source(ActionabilitySource.CIVIC)
                .build());

        evidenceItems.add(evidenceBuilder().event("BRAF p.Val600Glu")
                .drug("Sorafenib")
                .level(EvidenceLevel.LEVEL_B)
                .response("Responsive")
                .reference("variant:12")
                .source(ActionabilitySource.CIVIC)
                .build());

        evidenceItems.add(evidenceBuilder().event("BRAF p.Val600Glu")
                .drug("Trametinib + Dabrafenib")
                .level(EvidenceLevel.LEVEL_B)
                .response("Responsive")
                .reference("variant:12")
                .source(ActionabilitySource.CIVIC)
                .build());

        evidenceItems.add(evidenceBuilder().event("BRAF p.Val600Glu")
                .drug("Trametinib + Panitumumab + Dabrafenib")
                .level(EvidenceLevel.LEVEL_B)
                .response("Responsive")
                .reference("variant:12")
                .source(ActionabilitySource.CIVIC)
                .build());

        evidenceItems.add(evidenceBuilder().event("BRAF p.Val600Glu")
                .drug("Vemurafenib")
                .level(EvidenceLevel.LEVEL_B)
                .response("Resistant")
                .reference("variant:17")
                .source(ActionabilitySource.CIVIC)
                .build());

        evidenceItems.add(evidenceBuilder().event("BRAF p.Val600Glu")
                .drug("Vemurafenib + Cetuximab + Irinotecan")
                .level(EvidenceLevel.LEVEL_B)
                .response("Responsive")
                .reference("variant:12")
                .source(ActionabilitySource.CIVIC)
                .build());

        evidenceItems.add(evidenceBuilder().event("BRAF p.Val600Glu")
                .drug("Vemurafenib + Dabrafenib")
                .level(EvidenceLevel.LEVEL_B)
                .response("Responsive")
                .reference("variant:12")
                .source(ActionabilitySource.CIVIC)
                .build());

        evidenceItems.add(evidenceBuilder().event("PTEN Deletion")
                .drug("EGFR mAB inhibitor")
                .level(EvidenceLevel.LEVEL_B)
                .response("Resistant")
                .reference("PTEN:del")
                .source(ActionabilitySource.CGI)
                .build());

        evidenceItems.add(evidenceBuilder().event("PTEN Deletion")
                .drug("Everolimus")
                .level(EvidenceLevel.LEVEL_B)
                .response("Responsive")
                .reference("variant:213")
                .source(ActionabilitySource.CIVIC)
                .build());

        return evidenceItems;
    }

    @NotNull
    private static List<ClinicalTrial> createCOLO829ClinicalTrials() {
        List<ClinicalTrial> trials = Lists.newArrayList();
        trials.add(ImmutableClinicalTrial.builder()
                .event("BRAF p.Val600Glu")
                .acronym("IMPemBra")
                .reference("EXT8846 (NL54421.031.15)")
                .source(ActionabilitySource.ICLUSION)
                .isOnLabel(true)
                .build());

        trials.add(ImmutableClinicalTrial.builder()
                .event("BRAF p.Val600Glu")
                .acronym("PROCLAIM-001")
                .reference("EXT10151 (NL59299.042.17)")
                .source(ActionabilitySource.ICLUSION)
                .isOnLabel(true)
                .build());

        return trials;
    }

    @NotNull
    private static FittedPurity createFittedPurity(double impliedPurity, double averageTumorPloidy) {
        return ImmutableFittedPurity.builder()
                .purity(impliedPurity)
                .diploidProportion(0)
                .normFactor(0)
                .score(0)
                .ploidy(averageTumorPloidy)
                .somaticDeviation(0)
                .build();
    }

    @NotNull
    private static List<ReportableSomaticVariant> createCOLO829SomaticVariants(@NotNull PurityAdjuster purityAdjuster) {
        ReportableSomaticVariant variant1 = ImmutableReportableSomaticVariant.builder()
                .gene("BRAF")
                .isDrupActionable(true)
                .hgvsCodingImpact("c.1799T>A")
                .hgvsProteinImpact("p.Val600Glu")
                .hotspot(Hotspot.HOTSPOT)
                .clonality(Clonality.CLONAL)
                .alleleReadCount(107)
                .totalReadCount(161)
                .adjustedCopyNumber(6)
                .minorAllelePloidy(2)
                .biallelic(false)
                .adjustedVAF(purityAdjuster.purityAdjustedVAF("7", 6, 107D / 161D))
                .driverCategory(DriverCategory.ONCO)
                .driverLikelihood(1D)
                .build();

        ReportableSomaticVariant variant2 = ImmutableReportableSomaticVariant.builder()
                .gene("CDKN2A")
                .isDrupActionable(true)
                .hgvsCodingImpact("c.369_370delCG")
                .hgvsProteinImpact("p.Gly124fs")
                .hotspot(Hotspot.NEAR_HOTSPOT)
                .clonality(Clonality.CLONAL)
                .alleleReadCount(44)
                .totalReadCount(44)
                .adjustedCopyNumber(2)
                .minorAllelePloidy(0)
                .biallelic(true)
                .adjustedVAF(purityAdjuster.purityAdjustedVAF("9", 2, 44D / 44D))
                .driverCategory(DriverCategory.TSG)
                .driverLikelihood(0.9)
                .build();

        ReportableSomaticVariant variant3 = ImmutableReportableSomaticVariant.builder()
                .gene("SF3B1")
                .isDrupActionable(false)
                .hgvsCodingImpact("c.2153C>T")
                .hgvsProteinImpact("p.Pro718Leu")
                .hotspot(Hotspot.NON_HOTSPOT)
                .clonality(Clonality.CLONAL)
                .alleleReadCount(72)
                .totalReadCount(107)
                .adjustedCopyNumber(3)
                .minorAllelePloidy(1)
                .biallelic(false)
                .adjustedVAF(purityAdjuster.purityAdjustedVAF("2", 3, 72D / 107D))
                .driverCategory(DriverCategory.ONCO)
                .driverLikelihood(0.5)
                .build();

        ReportableSomaticVariant variant4 = ImmutableReportableSomaticVariant.builder()
                .gene("TP63")
                .isDrupActionable(false)
                .hgvsCodingImpact("c.1497G>T")
                .hgvsProteinImpact("p.Met499Ile")
                .hotspot(Hotspot.NON_HOTSPOT)
                .clonality(Clonality.CLONAL)
                .alleleReadCount(48)
                .totalReadCount(103)
                .adjustedCopyNumber(4)
                .minorAllelePloidy(2)
                .biallelic(false)
                .adjustedVAF(purityAdjuster.purityAdjustedVAF("3", 4, 48D / 103D))
                .driverCategory(DriverCategory.TSG)
                .driverLikelihood(0.1)
                .build();

        return Lists.newArrayList(variant1, variant2, variant3, variant4);
    }

    @NotNull
    private static List<GeneCopyNumber> createCOLO829CopyNumbers() {
        GeneCopyNumber copyNumber1 = createTestCopyNumberBuilder().chromosome("10")
                .chromosomeBand("q23.31")
                .gene("PTEN")
                .minCopyNumber(0)
                .maxCopyNumber(2)
                .build();

        return Lists.newArrayList(copyNumber1);
    }

    @NotNull
    private static List<ReportableGeneFusion> createCOLO829Fusions() {
        // KODU: Keep to be able to test actual fusions when needed.
        //        ReportableGeneFusion fusion1 = ImmutableReportableGeneFusion.builder()
        //                .geneStart("TMPRSS2")
        //                .geneStartTranscript("ENST00000398585")
        //                .geneContextStart("Intron 5")
        //                .geneEnd("PNPLA7")
        //                .geneEndTranscript("ENST00000406427")
        //                .geneContextEnd("Intron 3")
        //                .ploidy(0.4)
        //                .source(KnownFusionsModel.CIVIC)
        //                .build();
        //
        //        ReportableGeneFusion fusion2 = ImmutableReportableGeneFusion.builder()
        //                .geneStart("CLCN6")
        //                .geneStartTranscript("ENST00000346436")
        //                .geneContextStart("Intron 1")
        //                .geneEnd("BRAF")
        //                .geneEndTranscript("ENST00000288602")
        //                .geneContextEnd("Intron 8")
        //                .ploidy(1D)
        //                .source(KnownFusionsModel.ONCOKB)
        //                .build();
        //
        //        return Lists.newArrayList(fusion1, fusion2);
        return Lists.newArrayList();
    }

    @NotNull
    private static List<GermlineVariant> createCOLO829GermlineVariants(@NotNull PurityAdjuster purityAdjuster) {
        // KODU: Keep in case we need to generate germline variants for testing
        //        List<GermlineVariant> germlineVariants = Lists.newArrayList();
        //
        //        int totalReads = 112;
        //        int altReads = 67;
        //        double adjustedCopyNumber = 3D;
        //
        //        double adjustedVAF =
        //                purityAdjuster.purityAdjustedVAFWithHeterozygousNormal("13", adjustedCopyNumber, (double) altReads / (double) totalReads);
        //
        //        germlineVariants.add(ImmutableGermlineVariant.builder()
        //                .passFilter(true)
        //                .gene("BRCA2")
        //                .hgvsCodingImpact("c.5946delT")
        //                .hgvsProteinImpact("p.Ser1982fs")
        //                .totalReadCount(totalReads)
        //                .alleleReadCount(altReads)
        //                .germlineStatus("HET")
        //                .adjustedCopyNumber(adjustedCopyNumber)
        //                .adjustedVAF(adjustedVAF)
        //                .minorAllelePloidy(1D)
        //                .biallelic(false)
        //                .build());
        //
        //        return germlineVariants;
        return Lists.newArrayList();
    }

    @NotNull
    private static ChordAnalysis createCOLO829ChordAnalysis() {
        double brca1Value = 0D;
        double brca2Value = 0D;

        return ImmutableChordAnalysis.builder()
                .noneValue(1 - (brca1Value + brca2Value))
                .BRCA1Value(brca1Value)
                .BRCA2Value(brca2Value)
                .hrdValue(brca1Value + brca2Value)
                .predictedResponseValue(brca1Value + brca2Value > 0.5 ? 1 : 0)
                .build();
    }

    @NotNull
    private static List<ReportableGeneDisruption> createCOLO829Disruptions() {
        ReportableGeneDisruption disruption1 = createDisruptionBuilder().location("10q23.31")
                .gene("PTEN")
                .range("Intron 5 -> Intron 6")
                .type(StructuralVariantType.DEL)
                .ploidy(1.8)
                .geneMinCopies(0)
                .geneMaxCopies(2)
                .build();

        return Lists.newArrayList(disruption1);
    }

    @NotNull
    private static ImmutableEvidenceItem.Builder evidenceBuilder() {
        return ImmutableEvidenceItem.builder().drugsType(Strings.EMPTY).isOnLabel(false);
    }

    @NotNull
    private static ImmutableReportableGeneDisruption.Builder createDisruptionBuilder() {
        return ImmutableReportableGeneDisruption.builder().firstAffectedExon(1);
    }
}

