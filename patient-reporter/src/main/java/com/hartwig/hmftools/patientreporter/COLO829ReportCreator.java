package com.hartwig.hmftools.patientreporter;

import java.io.File;
import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;
import java.util.Locale;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.ActionabilitySource;
import com.hartwig.hmftools.common.actionability.ClinicalTrial;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.EvidenceLevel;
import com.hartwig.hmftools.common.actionability.EvidenceScope;
import com.hartwig.hmftools.common.actionability.ImmutableClinicalTrial;
import com.hartwig.hmftools.common.actionability.ImmutableEvidenceItem;
import com.hartwig.hmftools.common.center.Center;
import com.hartwig.hmftools.common.center.CenterModel;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ImmutableChordAnalysis;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.ecrf.projections.ImmutablePatientTumorLocation;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsFactory;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.ImmutableGeneCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.ImmutableFittedPurity;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.variant.Clonality;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.patientreporter.cfreport.CFReportWriter;
import com.hartwig.hmftools.patientreporter.germline.GermlineVariant;
import com.hartwig.hmftools.patientreporter.structural.ImmutableReportableGeneDisruption;
import com.hartwig.hmftools.patientreporter.structural.ReportableGeneDisruption;
import com.hartwig.hmftools.patientreporter.structural.ReportableGeneFusion;
import com.hartwig.hmftools.patientreporter.variants.ImmutableReportableSomaticVariant;
import com.hartwig.hmftools.patientreporter.variants.ReportableSomaticVariant;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class COLO829ReportCreator {

    private static final Logger LOGGER = LogManager.getLogger(PatientReporterApplication.class);

    private static final String SIGNATURE_PATH = "/signature/signature_test.png";
    private static final String RVA_LOGO_PATH = "/rva_logo/rva_logo_test.jpg";
    private static final String CENTER_CSV = "/csv/centers.csv";
    private static final String CENTER_MANUAL_CSV = "/csv/manual_mapping.csv";

    private static final String REPORT_BASE_DIR = System.getProperty("user.home") + "/hmf/tmp";

    private static final DateTimeFormatter DATE_FORMATTER = DateTimeFormatter.ofPattern("dd-MMM-yyyy", Locale.ENGLISH);

    public static void main(final String[] args) throws IOException {
        String resourceDir;

        if (args.length == 0) {
            resourceDir = System.getProperty("user.home") + "/hmf/repos/hmftools/patient-reporter/src/test/resources";
            LOGGER.info("Using default resource dir " + resourceDir);
        } else {
            resourceDir = args[0];
            LOGGER.info("Using configured resource dir " + resourceDir);
        }

        AnalysedPatientReport colo829Report = buildCOLO829(resourceDir);

        CFReportWriter reportWriter = new CFReportWriter();
        reportWriter.writeAnalysedPatientReport(colo829Report,
                getReportFilePath("hmf_test_sequence_report_" + String.valueOf(System.currentTimeMillis()) + ".pdf"));

    }

    @NotNull
    private static String getReportFilePath(@NotNull String filename) {
        return REPORT_BASE_DIR + File.separator + filename;
    }

    @NotNull
    private static AnalysedPatientReport buildCOLO829(@NotNull String resourceDir) {
        final double impliedTumorPurity = 1D;
        final double averageTumorPloidy = 3.1;
        final int tumorMutationalLoad = 182;
        final double tumorMutationalBurden = 13.6324;
        final double microsatelliteIndelsPerMb = 0.1179;

        final BaseReportData baseReportData = testBaseReportData(resourceDir);
        final FittedPurity fittedPurity = createFittedPurity(impliedTumorPurity, averageTumorPloidy);

        final PurityAdjuster purityAdjuster = new PurityAdjuster(Gender.MALE, fittedPurity);
        final List<EvidenceItem> tumorLocationSpecificEvidence = createCOLO829TumorSpecificEvidence();
        final List<ClinicalTrial> clinicalTrials = createCOLO829ClinicalTrials();
        final List<EvidenceItem> offLabelEvidence = createCOLO829OffLabelEvidence();
        final List<ReportableSomaticVariant> somaticVariants = createCOLO829SomaticVariants(purityAdjuster);
        final List<GeneCopyNumber> copyNumbers = createCOLO829CopyNumbers();
        final List<ReportableGeneFusion> fusions = Lists.newArrayList();
        final List<GermlineVariant> germlineVariants = Lists.newArrayList();
        final ChordAnalysis chordAnalysis = createCOLO829ChordAnalysis();
        final List<ReportableGeneDisruption> disruptions = createCOLO829Disruptions();

        final SampleReport sampleReport = createCOLO829SampleReport();

        return ImmutableAnalysedPatientReport.of(sampleReport,
                fittedPurity.purity(),
                true,
                fittedPurity.ploidy(),
                tumorLocationSpecificEvidence,
                clinicalTrials,
                offLabelEvidence,
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
                resourceDir + "/circos/circos_example.png",
                Optional.of("this is a test report and is based off COLO829"),
                baseReportData.signaturePath(),
                baseReportData.logoRVAPath());
    }

    @NotNull
    private static BaseReportData testBaseReportData(@NotNull String resourceDir) {
        try {
            final List<PatientTumorLocation> patientTumorLocations = Lists.newArrayList();
            final Lims lims = LimsFactory.empty();
            final CenterModel centerModel = Center.readFromCSV(resourceDir + CENTER_CSV, resourceDir + CENTER_MANUAL_CSV);
            return ImmutableBaseReportData.of(patientTumorLocations,
                    lims,
                    centerModel,
                    resourceDir + SIGNATURE_PATH,
                    resourceDir + RVA_LOGO_PATH);
        } catch (IOException exception) {
            throw new IllegalStateException("Could not generate test base reporter data: " + exception.getMessage());
        }
    }

    @NotNull
    private static SampleReport createCOLO829SampleReport() {
        final String sample = "COLO829T";
        return ImmutableSampleReport.of(sample,
                "A1",
                "A2",
                ImmutablePatientTumorLocation.of("COLO829", "Skin", "Melanoma"),
                Strings.EMPTY,
                "80%",
                LocalDate.parse("05-Jan-2018", DATE_FORMATTER),
                LocalDate.parse("01-Jan-2018", DATE_FORMATTER),
                "PREP013V23-QC037V20-SEQ008V25",
                "HMF Testing Center",
                "COLO",
                Strings.EMPTY,
                Strings.EMPTY,
                Strings.EMPTY,
                Strings.EMPTY);
    }

    @NotNull
    private static List<EvidenceItem> createCOLO829TumorSpecificEvidence() {
        List<EvidenceItem> evidenceItems = Lists.newArrayList();

        ImmutableEvidenceItem.Builder onLabelBuilder = evidenceBuilder().isOnLabel(true);

        evidenceItems.add(onLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Binimetinib + Encorafenib")
                .level(EvidenceLevel.LEVEL_A)
                .response("Responsive")
                .reference("V600E")
                .source(ActionabilitySource.ONCOKB)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(onLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Cobimetinib + Vemurafenib")
                .level(EvidenceLevel.LEVEL_A)
                .response("Responsive")
                .reference("V600E")
                .source(ActionabilitySource.ONCOKB)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(onLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Dabrafenib")
                .level(EvidenceLevel.LEVEL_A)
                .response("Responsive")
                .reference("V600E")
                .source(ActionabilitySource.ONCOKB)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(onLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Dabrafenib + Trametinib")
                .level(EvidenceLevel.LEVEL_A)
                .response("Responsive")
                .reference("V600E")
                .source(ActionabilitySource.ONCOKB)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(onLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Trametinib")
                .level(EvidenceLevel.LEVEL_A)
                .response("Responsive")
                .reference("V600E")
                .source(ActionabilitySource.ONCOKB)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(onLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Vemurafenib")
                .level(EvidenceLevel.LEVEL_A)
                .response("Responsive")
                .reference("V600E")
                .source(ActionabilitySource.ONCOKB)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(onLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Dabrafenib + Vemurafenib")
                .level(EvidenceLevel.LEVEL_B)
                .response("Responsive")
                .reference("variant:12")
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(onLabelBuilder.event("BRAF p.Val600Glu")
                .drug("RO4987655")
                .level(EvidenceLevel.LEVEL_B)
                .response("Responsive")
                .reference("variant:17")
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.BROAD)
                .build());

        return evidenceItems;
    }

    @NotNull
    private static List<ClinicalTrial> createCOLO829ClinicalTrials() {
        List<ClinicalTrial> trials = Lists.newArrayList();
        ImmutableClinicalTrial.Builder iclusionBuilder =
                ImmutableClinicalTrial.builder().cancerType(Strings.EMPTY).isOnLabel(true).source(ActionabilitySource.ICLUSION);

        trials.add(iclusionBuilder.event("BRAF p.Val600Glu")
                .scope(EvidenceScope.BROAD)
                .acronym("ERK inhibitor in tumors with MAPK pathway alterations")
                .reference("EXT10454 (NL57739.031.16)")
                .build());
        trials.add(iclusionBuilder.event("BRAF p.Val600Glu")
                .scope(EvidenceScope.SPECIFIC)
                .acronym("IMPemBra")
                .reference("EXT8846 (NL54421.031.15)")
                .build());
        trials.add(iclusionBuilder.event("BRAF p.Val600Glu")
                .scope(EvidenceScope.BROAD)
                .acronym("LXH254 in tumors with MAPK pathway alterations")
                .reference("EXT10453 (NL55506.078.15)")
                .build());
        trials.add(iclusionBuilder.event("BRAF p.Val600Glu")
                .scope(EvidenceScope.BROAD)
                .acronym("Novartis CTMT212X2102")
                .reference("EXT3437 (NL56240.056.16)")
                .build());
        trials.add(iclusionBuilder.event("BRAF p.Val600Glu")
                .scope(EvidenceScope.SPECIFIC)
                .acronym("PROCLAIM-001")
                .reference("EXT10151 (NL59299.042.17)")
                .build());
        trials.add(iclusionBuilder.event("BRAF p.Val600Glu")
                .scope(EvidenceScope.BROAD)
                .acronym("iMATRIXcobi (GO29665)")
                .reference("EXT6769 (NL52503.078.16)")
                .build());
        trials.add(iclusionBuilder.event("PTEN Deletion")
                .scope(EvidenceScope.SPECIFIC)
                .acronym("AZD5363 (AKT inhibitor) in advanced solid malignancies")
                .reference("EXT2552 (NL33755.031.10)")
                .build());

        return trials;
    }

    @NotNull
    private static List<EvidenceItem> createCOLO829OffLabelEvidence() {
        List<EvidenceItem> evidenceItems = Lists.newArrayList();

        ImmutableEvidenceItem.Builder offLabelBuilder = evidenceBuilder().isOnLabel(false);

        evidenceItems.add(offLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Alpelisib + Cetuximab + Encorafenib")
                .level(EvidenceLevel.LEVEL_B)
                .response("Responsive")
                .reference("variant:17")
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.BROAD)
                .build());

        evidenceItems.add(offLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Bevacizumab")
                .level(EvidenceLevel.LEVEL_B)
                .response("Resistant")
                .reference("variant:12")
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(offLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Cetuximab")
                .level(EvidenceLevel.LEVEL_B)
                .response("Resistant")
                .reference("BRAF:V600E")
                .source(ActionabilitySource.CGI)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(offLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Cetuximab + Encorafenib")
                .level(EvidenceLevel.LEVEL_B)
                .response("Responsive")
                .reference("variant:17")
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.BROAD)
                .build());

        evidenceItems.add(offLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Cetuximab + Irinotecan + Vemurafenib")
                .level(EvidenceLevel.LEVEL_B)
                .response("Responsive")
                .reference("variant:12")
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(offLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Dabrafenib + Panitumumab + Trametinib")
                .level(EvidenceLevel.LEVEL_B)
                .response("Responsive")
                .reference("variant:12")
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(offLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Irinotecan")
                .level(EvidenceLevel.LEVEL_B)
                .response("Resistant")
                .reference("variant:399")
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(offLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Oxaliplatin")
                .level(EvidenceLevel.LEVEL_B)
                .response("Resistant")
                .reference("variant:12")
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(offLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Panitumumab")
                .level(EvidenceLevel.LEVEL_B)
                .response("Resistant")
                .reference("BRAF:V600E")
                .source(ActionabilitySource.CGI)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(offLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Sorafenib")
                .level(EvidenceLevel.LEVEL_B)
                .response("Responsive")
                .reference("variant:12")
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(offLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Vemurafenib")
                .level(EvidenceLevel.LEVEL_B)
                .response("Resistant")
                .reference("variant:17")
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.BROAD)
                .build());

        evidenceItems.add(offLabelBuilder.event("PTEN Deletion")
                .drug("EGFR mAB inhibitor")
                .level(EvidenceLevel.LEVEL_B)
                .response("Resistant")
                .reference("PTEN:del")
                .source(ActionabilitySource.CGI)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(offLabelBuilder.event("PTEN Deletion")
                .drug("Everolimus")
                .level(EvidenceLevel.LEVEL_B)
                .response("Responsive")
                .reference("variant:213")
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        return evidenceItems;
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
    private static ImmutableGeneCopyNumber.Builder createTestCopyNumberBuilder() {
        return ImmutableGeneCopyNumber.builder()
                .start(1)
                .end(2)
                .gene(Strings.EMPTY)
                .chromosome("1")
                .chromosomeBand(Strings.EMPTY)
                .minRegionStart(0)
                .minRegionStartSupport(SegmentSupport.NONE)
                .minRegionEnd(0)
                .minRegionEndSupport(SegmentSupport.NONE)
                .minRegionMethod(CopyNumberMethod.UNKNOWN)
                .minRegions(1)
                .germlineHet2HomRegions(0)
                .germlineHomRegions(0)
                .somaticRegions(1)
                .minCopyNumber(0)
                .maxCopyNumber(0)
                .transcriptID("trans")
                .transcriptVersion(0)
                .nonsenseBiallelicCount(0)
                .nonsenseNonBiallelicCount(0)
                .nonsenseNonBiallelicPloidy(0)
                .spliceBiallelicCount(0)
                .spliceNonBiallelicCount(0)
                .spliceNonBiallelicPloidy(0)
                .missenseBiallelicCount(0)
                .missenseNonBiallelicCount(0)
                .missenseNonBiallelicPloidy(0)
                .minMinorAllelePloidy(0);
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
                .type("DEL")
                .ploidy(2D)
                .geneMinCopies(0)
                .geneMaxCopies(2)
                .build();

        return Lists.newArrayList(disruption1);
    }

    @NotNull
    private static ImmutableEvidenceItem.Builder evidenceBuilder() {
        return ImmutableEvidenceItem.builder().drugsType(Strings.EMPTY).cancerType(Strings.EMPTY);
    }

    @NotNull
    private static ImmutableReportableGeneDisruption.Builder createDisruptionBuilder() {
        return ImmutableReportableGeneDisruption.builder().firstAffectedExon(1);
    }
}
