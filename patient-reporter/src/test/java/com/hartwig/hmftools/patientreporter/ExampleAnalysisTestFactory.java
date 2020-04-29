package com.hartwig.hmftools.patientreporter;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testReportData;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;
import java.util.Locale;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.actionability.ActionabilitySource;
import com.hartwig.hmftools.common.actionability.ClinicalTrial;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.EvidenceLevel;
import com.hartwig.hmftools.common.actionability.EvidenceScope;
import com.hartwig.hmftools.common.actionability.ImmutableClinicalTrial;
import com.hartwig.hmftools.common.actionability.ImmutableEvidenceItem;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ImmutableChordAnalysis;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.ecrf.projections.ImmutablePatientTumorLocation;
import com.hartwig.hmftools.common.lims.LimsStudy;
import com.hartwig.hmftools.common.lims.LimsViralInsertionChoice;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.structural.annotation.ImmutableReportableGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusion;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;
import com.hartwig.hmftools.patientreporter.copynumber.ImmutableReportableGainLoss;
import com.hartwig.hmftools.patientreporter.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.patientreporter.homozygousdisruption.ImmutableReportableHomozygousDisruption;
import com.hartwig.hmftools.patientreporter.homozygousdisruption.ReportableHomozygousDisruption;
import com.hartwig.hmftools.patientreporter.qcfail.ImmutableQCFailReport;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReport;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailStudy;
import com.hartwig.hmftools.patientreporter.structural.ImmutableReportableGeneDisruption;
import com.hartwig.hmftools.patientreporter.structural.ReportableGeneDisruption;
import com.hartwig.hmftools.patientreporter.variants.ImmutableReportableVariant;
import com.hartwig.hmftools.patientreporter.variants.ReportableVariant;
import com.hartwig.hmftools.patientreporter.viralInsertion.ImmutableViralInsertion;
import com.hartwig.hmftools.patientreporter.viralInsertion.ViralInsertion;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ExampleAnalysisTestFactory {

    private static final DateTimeFormatter DATE_FORMATTER = DateTimeFormatter.ofPattern("dd-MMM-yyyy", Locale.ENGLISH);
    private static final String CIRCOS_PATH = Resources.getResource("test_run/purple/plot/sample.circos.png").getPath();

    private ExampleAnalysisTestFactory() {
    }

    @NotNull
    public static AnalysedPatientReport buildCOLO829() {
        final boolean hasReliablePurity = true;
        final double impliedTumorPurity = 1D;
        final double averageTumorPloidy = 3.1;
        final int tumorMutationalLoad = 180;
        final double tumorMutationalBurden = 13.6;
        final double microsatelliteIndelsPerMb = 0.11;
        final MicrosatelliteStatus msiStatus = MicrosatelliteStatus.MSS;
        final TumorMutationalStatus tmlStatus = TumorMutationalStatus.HIGH;

        final ReportData reportData = testReportData();

        final List<EvidenceItem> tumorLocationSpecificEvidence = createCOLO829TumorSpecificEvidence();
        final List<ClinicalTrial> clinicalTrials = createCOLO829ClinicalTrials();
        final List<EvidenceItem> offLabelEvidence = createCOLO829OffLabelEvidence();
        final List<ReportableVariant> reportableVariants = createCOLO829SomaticVariants();
        final List<ReportableGainLoss> gainsAndLosses = createCOLO829GainsLosses();
        final List<ReportableGeneFusion> fusions = Lists.newArrayList();
        final List<ReportableHomozygousDisruption> reportableHomozygousDisruptions = Lists.newArrayList();
        final List<ReportableGeneDisruption> disruptions = createCOLO829Disruptions();
        final ChordAnalysis chordAnalysis = createCOLO829ChordAnalysis();
        final List<ViralInsertion> viralInsertions = Lists.newArrayList();

        final String sampleId = "PNT00012345T";
        final SampleReport sampleReport = createSkinMelanomaSampleReport(sampleId);

        final String clinicalSummary = "Melanoma sample showing:\n"
                + " - activating BRAF mutation that is associated with response to BRAF-inhibitors (in combination with a MEK-inhibitor)\n"
                + " - complete inactivation of CDKN2A, indicating potential benefit of CDK4/6 inhibitors\n"
                + " - complete inactivation/loss of PTEN likely resulting in an activation of the PI3K-AKT-mTOR pathway "
                + "and indicating potential benefit of mTOR/PI3K inhibitors\n"
                + " - high mutational burden (mutational load (ML) of 180, tumor mutation burden (TMB) of 13.6) that is "
                + "potentially associated with an increased response rate to checkpoint inhibitor immunotherapy";

        return ImmutableAnalysedPatientReport.builder()
                .sampleReport(sampleReport)
                .impliedPurity(impliedTumorPurity)
                .hasReliablePurity(hasReliablePurity)
                .hasReliableQuality(true)
                .averageTumorPloidy(averageTumorPloidy)
                .reportableViralInsertions(LimsViralInsertionChoice.NO_REPORT_VIRAL_INSERTIONS)
                .clinicalSummary(clinicalSummary)
                .tumorSpecificEvidence(tumorLocationSpecificEvidence)
                .clinicalTrials(clinicalTrials)
                .offLabelEvidence(offLabelEvidence)
                .reportableVariants(reportableVariants)
                .microsatelliteIndelsPerMb(microsatelliteIndelsPerMb)
                .microsatelliteStatus(msiStatus)
                .tumorMutationalLoad(tumorMutationalLoad)
                .tumorMutationalLoadStatus(tmlStatus)
                .tumorMutationalBurden(tumorMutationalBurden)
                .chordAnalysis(chordAnalysis)
                .gainsAndLosses(gainsAndLosses)
                .geneFusions(fusions)
                .geneDisruptions(disruptions)
                .reportableHomozygousDisruptions(reportableHomozygousDisruptions)
                .viralInsertions(viralInsertions)
                .circosPath(CIRCOS_PATH)
                .comments(Optional.of("This is a test report and is based off COLO829"))
                .isCorrectedReport(false)
                .isUnofficialReport(false)
                .signaturePath(reportData.signaturePath())
                .logoRVAPath(reportData.logoRVAPath())
                .logoCompanyPath(reportData.logoCompanyPath())
                .build();
    }

    @NotNull
    public static AnalysedPatientReport buildAnalysisWithAllTablesFilledIn(@NotNull String sampleId) {
        final boolean hasReliablePurity = true;
        final double impliedTumorPurity = 1D;
        final double averageTumorPloidy = 3.1;
        final int tumorMutationalLoad = 182;
        final double tumorMutationalBurden = 13.6;
        final double microsatelliteIndelsPerMb = 0.1089;
        final MicrosatelliteStatus msiStatus = MicrosatelliteStatus.MSS;
        final TumorMutationalStatus tmlStatus = TumorMutationalStatus.HIGH;

        final ReportData reportData = testReportData();

        final List<EvidenceItem> tumorLocationSpecificEvidence = createCOLO829TumorSpecificEvidence();
        final List<ClinicalTrial> clinicalTrials = createCOLO829ClinicalTrials();
        final List<EvidenceItem> offLabelEvidence = createCOLO829OffLabelEvidence();
        final List<ReportableVariant> reportableVariants = createAllSomaticVariants();
        final List<ReportableGainLoss> gainsAndLosses = createCOLO829GainsLosses();
        final List<ReportableGeneFusion> fusions = createTestFusions();
        final ChordAnalysis chordAnalysis = createCOLO829ChordAnalysis();
        final List<ReportableGeneDisruption> disruptions = createCOLO829Disruptions();
        final List<ViralInsertion> viralInsertions = createTestViralInsertions();
        final List<ReportableHomozygousDisruption> reportableHomozygousDisruptions = createTestHomozygousDisruptions();

        final SampleReport sampleReport = createSkinMelanomaSampleReport(sampleId);
        final String clinicalSummary = Strings.EMPTY;

        return ImmutableAnalysedPatientReport.builder()
                .sampleReport(sampleReport)
                .impliedPurity(impliedTumorPurity)
                .hasReliablePurity(hasReliablePurity)
                .hasReliableQuality(true)
                .averageTumorPloidy(averageTumorPloidy)
                .reportableViralInsertions(LimsViralInsertionChoice.REPORT_VIRAL_INSERION)
                .clinicalSummary(clinicalSummary)
                .tumorSpecificEvidence(tumorLocationSpecificEvidence)
                .clinicalTrials(clinicalTrials)
                .offLabelEvidence(offLabelEvidence)
                .reportableVariants(reportableVariants)
                .microsatelliteIndelsPerMb(microsatelliteIndelsPerMb)
                .microsatelliteStatus(msiStatus)
                .tumorMutationalLoad(tumorMutationalLoad)
                .tumorMutationalLoadStatus(tmlStatus)
                .tumorMutationalBurden(tumorMutationalBurden)
                .chordAnalysis(chordAnalysis)
                .gainsAndLosses(gainsAndLosses)
                .geneFusions(fusions)
                .geneDisruptions(disruptions)
                .reportableHomozygousDisruptions(reportableHomozygousDisruptions)
                .viralInsertions(viralInsertions)
                .circosPath(CIRCOS_PATH)
                .comments(Optional.of("This is a test report and does not relate to any real patient"))
                .isCorrectedReport(false)
                .isUnofficialReport(false)
                .signaturePath(reportData.signaturePath())
                .logoRVAPath(reportData.logoRVAPath())
                .logoCompanyPath(reportData.logoCompanyPath())
                .build();
    }

    @NotNull
    public static AnalysedPatientReport buildAnalysisWithAllTablesForBelowDetectionLimitSample(@NotNull String sampleId) {
        final boolean hasReliablePurity = false;
        final double impliedTumorPurity = 1D;
        final double averageTumorPloidy = 3.1;
        final int tumorMutationalLoad = 182;
        final double tumorMutationalBurden = 13.6;
        final double microsatelliteIndelsPerMb = 0.1089;
        final MicrosatelliteStatus msiStatus = MicrosatelliteStatus.MSS;
        final TumorMutationalStatus tmlStatus = TumorMutationalStatus.HIGH;

        final ReportData reportData = testReportData();

        final List<EvidenceItem> tumorLocationSpecificEvidence = createCOLO829TumorSpecificEvidence();
        final List<ClinicalTrial> clinicalTrials = createCOLO829ClinicalTrials();
        final List<EvidenceItem> offLabelEvidence = createCOLO829OffLabelEvidence();
        final List<ReportableVariant> reportableVariants = createAllSomaticVariants();
        final List<ReportableGainLoss> gainsAndLosses = createCOLO829GainsLosses();
        final List<ReportableGeneFusion> fusions = createTestFusions();
        final ChordAnalysis chordAnalysis = createCOLO829ChordAnalysis();
        final List<ReportableGeneDisruption> disruptions = createCOLO829Disruptions();
        final List<ViralInsertion> viralInsertions = createTestViralInsertions();
        final List<ReportableHomozygousDisruption> reportableHomozygousDisruptions = createTestHomozygousDisruptions();

        final SampleReport sampleReport = createSkinMelanomaSampleReport(sampleId);
        final String clinicalSummary = Strings.EMPTY;

        return ImmutableAnalysedPatientReport.builder()
                .sampleReport(sampleReport)
                .impliedPurity(impliedTumorPurity)
                .hasReliablePurity(hasReliablePurity)
                .hasReliableQuality(true)
                .averageTumorPloidy(averageTumorPloidy)
                .reportableViralInsertions(LimsViralInsertionChoice.REPORT_VIRAL_INSERION)
                .clinicalSummary(clinicalSummary)
                .tumorSpecificEvidence(tumorLocationSpecificEvidence)
                .clinicalTrials(clinicalTrials)
                .offLabelEvidence(offLabelEvidence)
                .reportableVariants(reportableVariants)
                .microsatelliteIndelsPerMb(microsatelliteIndelsPerMb)
                .microsatelliteStatus(msiStatus)
                .tumorMutationalLoad(tumorMutationalLoad)
                .tumorMutationalLoadStatus(tmlStatus)
                .tumorMutationalBurden(tumorMutationalBurden)
                .chordAnalysis(chordAnalysis)
                .gainsAndLosses(gainsAndLosses)
                .geneFusions(fusions)
                .geneDisruptions(disruptions)
                .reportableHomozygousDisruptions(reportableHomozygousDisruptions)
                .viralInsertions(viralInsertions)
                .circosPath(CIRCOS_PATH)
                .comments(Optional.of("This is a test report and does not relate to any real patient"))
                .isCorrectedReport(false)
                .isUnofficialReport(false)
                .signaturePath(reportData.signaturePath())
                .logoRVAPath(reportData.logoRVAPath())
                .logoCompanyPath(reportData.logoCompanyPath())
                .build();
    }

    public static QCFailReport buildQCFailReport(@NotNull String sampleId, @NotNull QCFailReason reason) {
        SampleReport sampleReport = createSkinMelanomaSampleReport(sampleId);

        LimsStudy study = LimsStudy.fromSampleId(sampleId);
        QCFailStudy failStudy;
        switch (study) {
            case CORE:
                failStudy = QCFailStudy.CORE;
                break;
            case WIDE:
                failStudy = QCFailStudy.WIDE;
                break;
            case DRUP:
                failStudy = QCFailStudy.DRUP;
                break;
            default:
                failStudy = QCFailStudy.CPCT;
        }

        final ReportData reportData = testReportData();
        return ImmutableQCFailReport.builder()
                .sampleReport(sampleReport)
                .reason(reason)
                .study(failStudy)
                .comments(Optional.empty())
                .isCorrectedReport(false)
                .signaturePath(reportData.signaturePath())
                .logoRVAPath(reportData.logoRVAPath())
                .logoCompanyPath(reportData.logoCompanyPath())
                .build();
    }

    @NotNull
    private static SampleReport createSkinMelanomaSampleReport(@NotNull String sample) {
        SampleMetadata sampleMetadata = ImmutableSampleMetadata.builder()
                .refSampleId(Strings.EMPTY)
                .refSampleBarcode("FR12123488")
                .tumorSampleId(sample)
                .tumorSampleBarcode("FR12345678")
                .build();

        return ImmutableSampleReport.builder()
                .sampleMetadata(sampleMetadata)
                .patientTumorLocation(ImmutablePatientTumorLocation.of(Strings.EMPTY, "Skin", "Melanoma"))
                .refArrivalDate(LocalDate.parse("01-Jan-2019", DATE_FORMATTER))
                .tumorArrivalDate(LocalDate.parse("05-Jan-2019", DATE_FORMATTER))
                .purityShallowSeq(Strings.EMPTY)
                .labProcedures("PREP013V23-QC037V20-SEQ008V25")
                .requesterName("C")
                .requesterEmail("D")
                .addressee("HMF Testing Center")
                .hospitalName(Strings.EMPTY)
                .hospitalPIName(Strings.EMPTY)
                .hospitalPIEmail(Strings.EMPTY)
                .projectName("TEST")
                .submissionId("10")
                .hospitalPatientId("4567")
                .hospitalPathologySampleId("1234")
                .studyRequesterName("contact")
                .studyRequesterEmail("contact@.nl")
                .build();
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
                .drug("RO4987655")
                .level(EvidenceLevel.LEVEL_B)
                .response("Responsive")
                .reference("variant:208")
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.GENE_LEVEL)
                .build());

        return evidenceItems;
    }

    @NotNull
    private static List<ClinicalTrial> createCOLO829ClinicalTrials() {
        List<ClinicalTrial> trials = Lists.newArrayList();
        ImmutableClinicalTrial.Builder iclusionBuilder =
                ImmutableClinicalTrial.builder().cancerType(Strings.EMPTY).isOnLabel(true).source(ActionabilitySource.ICLUSION);

        trials.add(iclusionBuilder.event("BRAF p.Val600Glu")
                .scope(EvidenceScope.GENE_LEVEL)
                .acronym("DRUP")
                .reference("EXT10299 (NL54757.031.16)")
                .build());
        trials.add(iclusionBuilder.event("BRAF p.Val600Glu")
                .scope(EvidenceScope.GENE_LEVEL)
                .acronym("EBIN (EORTC-1612-MG)")
                .reference("EXT11284 (NL67202.031.18)")
                .build());
        trials.add(iclusionBuilder.event("BRAF p.Val600Glu")
                .scope(EvidenceScope.GENE_LEVEL)
                .acronym("LXH254 in tumors with MAPK pathway alterations")
                .reference("EXT10453 (NL55506.078.15)")
                .build());
        trials.add(iclusionBuilder.event("BRAF p.Val600Glu")
                .scope(EvidenceScope.GENE_LEVEL)
                .acronym("POLARIS")
                .reference("EXT11388 (NL69569.028.19)")
                .build());
        trials.add(iclusionBuilder.event("BRAF p.Val600Glu")
                .scope(EvidenceScope.SPECIFIC)
                .acronym("PROCLAIM-001")
                .reference("EXT10241 (NL59299.042.17)")
                .build());
        trials.add(iclusionBuilder.event("BRAF p.Val600Glu")
                .scope(EvidenceScope.SPECIFIC)
                .acronym("REDUCTOR")
                .reference("EXT6690 (NL45261.031.13)")
                .build());
        trials.add(iclusionBuilder.event("CDKN2A p.Ala68fs")
                .scope(EvidenceScope.GENE_LEVEL)
                .acronym("DRUP")
                .reference("EXT10299 (NL54757.031.16)")
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
                .scope(EvidenceScope.GENE_LEVEL)
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
                .drug("CI-1040")
                .level(EvidenceLevel.LEVEL_B)
                .response("Responsive")
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
                .scope(EvidenceScope.GENE_LEVEL)
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
                .reference("variant:12")
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
                .drug("Vemurafenib")
                .level(EvidenceLevel.LEVEL_B)
                .response("Resistant")
                .reference("variant:17")
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.GENE_LEVEL)
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
    private static List<ReportableVariant> createCOLO829SomaticVariants() {
        ReportableVariant variant1 = ImmutableReportableVariant.builder()
                .gene("BRAF")
                .position(140453136)
                .chromosome("7")
                .ref("T")
                .alt("A")
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .canonicalHgvsCodingImpact("c.1799T>A")
                .canonicalHgvsProteinImpact("p.Val600Glu")
                .notifyClinicalGeneticist(false)
                .driverCategory(DriverCategory.ONCO)
                .gDNA("7:140453136")
                .alleleReadCount(154)
                .totalReadCount(225)
                .allelePloidy(4.1)
                .totalPloidy(6.1)
                .hotspot(Hotspot.HOTSPOT)
                .biallelic(false)
                .driverLikelihood(1D)
                .clonalLikelihood(1D)
                .build();

        ReportableVariant variant2 = ImmutableReportableVariant.builder()
                .gene("CDKN2A")
                .position(21971153)
                .chromosome("9")
                .ref("CCG")
                .alt("C")
                .canonicalCodingEffect(CodingEffect.NONSENSE_OR_FRAMESHIFT)
                .canonicalHgvsCodingImpact("c.203_204delCG")
                .canonicalHgvsProteinImpact("p.Ala68fs")
                .notifyClinicalGeneticist(false)
                .driverCategory(DriverCategory.TSG)
                .gDNA("9:21971153")
                .alleleReadCount(95)
                .totalReadCount(95)
                .allelePloidy(2D)
                .totalPloidy(2D)
                .hotspot(Hotspot.NEAR_HOTSPOT)
                .biallelic(true)
                .clonalLikelihood(1D)
                .driverLikelihood(0.9)
                .build();

        ReportableVariant variant3 = ImmutableReportableVariant.builder()
                .gene("TERT")
                .position(1295228)
                .chromosome("5")
                .ref("GG")
                .alt("AA")
                .canonicalCodingEffect(CodingEffect.NONE)
                .canonicalHgvsCodingImpact("c.-125_-124delCCinsTT")
                .canonicalHgvsProteinImpact(Strings.EMPTY)
                .notifyClinicalGeneticist(false)
                .driverCategory(DriverCategory.ONCO)
                .gDNA("5:1295228")
                .alleleReadCount(49)
                .totalReadCount(49)
                .allelePloidy(2D)
                .totalPloidy(2D)
                .hotspot(Hotspot.HOTSPOT)
                .biallelic(false)
                .clonalLikelihood(1D)
                .driverLikelihood(0.85)
                .build();

        ReportableVariant variant4 = ImmutableReportableVariant.builder()
                .gene("SF3B1")
                .position(198266779)
                .chromosome("2")
                .ref("C")
                .alt("T")
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .canonicalHgvsCodingImpact("c.2153C>T")
                .canonicalHgvsProteinImpact("p.Pro718Leu")
                .notifyClinicalGeneticist(false)
                .driverCategory(DriverCategory.ONCO)
                .gDNA("2:198266779")
                .alleleReadCount(76)
                .totalReadCount(115)
                .allelePloidy(2D)
                .totalPloidy(3.1)
                .hotspot(Hotspot.NON_HOTSPOT)
                .biallelic(false)
                .clonalLikelihood(1D)
                .driverLikelihood(0.5)
                .build();

        ReportableVariant variant5 = ImmutableReportableVariant.builder()
                .gene("TP63")
                .position(189604330)
                .chromosome("3")
                .ref("G")
                .alt("T")
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .canonicalHgvsCodingImpact("c.1497G>T")
                .canonicalHgvsProteinImpact("p.Met499Ile")
                .notifyClinicalGeneticist(false)
                .driverCategory(DriverCategory.TSG)
                .gDNA("3:189604330")
                .alleleReadCount(52)
                .totalReadCount(119)
                .allelePloidy(1.8)
                .totalPloidy(4D)
                .hotspot(Hotspot.NON_HOTSPOT)
                .biallelic(false)
                .clonalLikelihood(1D)
                .driverLikelihood(0.1)
                .build();

        return Lists.newArrayList(variant1, variant2, variant3, variant4, variant5);
    }

    @NotNull
    private static List<ReportableVariant> createAllSomaticVariants() {
        ReportableVariant variant1 = ImmutableReportableVariant.builder()
                .gene("TP63")
                .position(189604330)
                .chromosome("3")
                .ref("G")
                .alt("T")
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .canonicalHgvsCodingImpact("c.1497G>T")
                .canonicalHgvsProteinImpact("p.Met499Ile")
                .notifyClinicalGeneticist(false)
                .driverCategory(DriverCategory.TSG)
                .gDNA("3:189604330")
                .alleleReadCount(48)
                .totalReadCount(103)
                .allelePloidy(2.1)
                .totalPloidy(4.1)
                .biallelic(false)
                .hotspot(Hotspot.NON_HOTSPOT)
                .clonalLikelihood(0.47)
                .driverLikelihood(0.1)
                .build();

        ReportableVariant variant2 = ImmutableReportableVariant.builder()
                .gene("KIT")
                .position(81627197)
                .chromosome("3")
                .ref("G")
                .alt("T")
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .canonicalHgvsCodingImpact("c.1497G>T")
                .canonicalHgvsProteinImpact("p.Met499Ile")
                .notifyClinicalGeneticist(true)
                .driverCategory(DriverCategory.TSG)
                .gDNA("3:81627197")
                .alleleReadCount(48)
                .totalReadCount(103)
                .allelePloidy(1.3)
                .totalPloidy(2.5)
                .hotspot(Hotspot.NON_HOTSPOT)
                .biallelic(true)
                .clonalLikelihood(0.68)
                .driverLikelihood(0.1)
                .build();

        return Lists.newArrayList(variant1, variant2);
    }

    @NotNull
    private static List<ReportableGainLoss> createCOLO829GainsLosses() {
        ReportableGainLoss gainLoss1 = ImmutableReportableGainLoss.builder()
                .chromosome("10")
                .chromosomeBand("q23.31")
                .gene("PTEN")
                .copies(0)
                .interpretation(CopyNumberInterpretation.PARTIAL_LOSS)
                .build();

        return Lists.newArrayList(gainLoss1);
    }

    @NotNull
    private static List<ReportableGeneFusion> createTestFusions() {
        ReportableGeneFusion fusion1 = ImmutableReportableGeneFusion.builder()
                .geneStart("TMPRSS2")
                .geneTranscriptStart("ENST00000398585")
                .geneContextStart("Intron 5")
                .geneEnd("PNPLA7")
                .geneTranscriptEnd("ENST00000406427")
                .geneContextEnd("Intron 3")
                .ploidy(0.4)
                .build();

        ReportableGeneFusion fusion2 = ImmutableReportableGeneFusion.builder()
                .geneStart("CLCN6")
                .geneTranscriptStart("ENST00000346436")
                .geneContextStart("Intron 1")
                .geneEnd("BRAF")
                .geneTranscriptEnd("ENST00000288602")
                .geneContextEnd("Intron 8")
                .ploidy(1D)
                .build();

        return Lists.newArrayList(fusion1, fusion2);
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
                .predictedResponseValue(brca1Value + brca2Value > 0.5)
                .build();
    }

    @NotNull
    private static List<ReportableGeneDisruption> createCOLO829Disruptions() {
        ReportableGeneDisruption disruption1 = createDisruptionBuilder().location("10q23.31")
                .gene("PTEN")
                .range("Intron 5 -> Intron 6")
                .type("DEL")
                .ploidy(2D)
                .undisruptedCopyNumber(0)
                .build();

        return Lists.newArrayList(disruption1);
    }

    @NotNull
    private static List<ReportableHomozygousDisruption> createTestHomozygousDisruptions() {
        List<ReportableHomozygousDisruption> homozygousDisruptions = Lists.newArrayList(ImmutableReportableHomozygousDisruption.builder()
                .chromosome("8")
                .chromosomeBand("p22")
                .gene("SGCZ")
                .build());
        return Lists.newArrayList(homozygousDisruptions);
    }

    @NotNull
    private static List<ViralInsertion> createTestViralInsertions() {
        List<ViralInsertion> viralInsertions =
                Lists.newArrayList(ImmutableViralInsertion.builder().virus("Human papillomavirus type 16").viralInsertionCount(2).build());
        return Lists.newArrayList(viralInsertions);
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

