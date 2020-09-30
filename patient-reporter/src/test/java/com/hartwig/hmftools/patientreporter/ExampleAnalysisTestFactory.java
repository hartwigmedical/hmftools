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
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.clinical.ImmutablePatientTumorLocation;
import com.hartwig.hmftools.common.fusion.ImmutableReportableGeneFusion;
import com.hartwig.hmftools.common.fusion.ReportableGeneFusion;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.hospital.HospitalContactData;
import com.hartwig.hmftools.common.lims.hospital.ImmutableHospitalContactData;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.copynumber.ImmutableReportableGainLoss;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;
import com.hartwig.hmftools.patientreporter.homozygousdisruption.ImmutableReportableHomozygousDisruption;
import com.hartwig.hmftools.patientreporter.homozygousdisruption.ReportableHomozygousDisruption;
import com.hartwig.hmftools.patientreporter.qcfail.ImmutableQCFailReport;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReport;
import com.hartwig.hmftools.patientreporter.structural.ImmutableReportableGeneDisruption;
import com.hartwig.hmftools.patientreporter.structural.ReportableGeneDisruption;
import com.hartwig.hmftools.patientreporter.variants.ImmutableReportableVariant;
import com.hartwig.hmftools.patientreporter.variants.ReportableVariant;
import com.hartwig.hmftools.patientreporter.viralInsertion.ImmutableViralInsertion;
import com.hartwig.hmftools.patientreporter.viralInsertion.ViralInsertion;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ExampleAnalysisTestFactory {

    private static final DateTimeFormatter DATE_FORMATTER = DateTimeFormatter.ofPattern("dd-MMM-yyyy", Locale.ENGLISH);
    private static final String CIRCOS_PATH = Resources.getResource("test_run/purple/plot/sample.circos.png").getPath();

    private ExampleAnalysisTestFactory() {
    }

    @NotNull
    public static AnalysedPatientReport buildCOLO829(@NotNull String sampleId, boolean correctionReport, @Nullable String comments) {
        return buildWithCOLO829Data(sampleId, correctionReport, comments, QsFormNumber.FOR_080.display(), true, 1D, true, false);
    }

    @NotNull
    public static AnalysedPatientReport buildWithCOLO829Data(@NotNull String sampleId, boolean correctionReport, @Nullable String comments,
            @NotNull String qcForNumber, boolean hasReliablePurity, double impliedTumorPurity, boolean includeSummary, boolean reportGermline) {
        double averageTumorPloidy = 3.1;
        int tumorMutationalLoad = 190;
        double tumorMutationalBurden = 13.7;
        double microsatelliteIndelsPerMb = 0.12;
        double chordHrdValue = 0D;

        ReportData reportData = testReportData();

        List<EvidenceItem> tumorLocationSpecificEvidence = createCOLO829TumorSpecificEvidence();
        List<ClinicalTrial> clinicalTrials = createCOLO829ClinicalTrials();
        List<EvidenceItem> offLabelEvidence = createCOLO829OffLabelEvidence();
        List<ReportableVariant> reportableVariants = createCOLO829SomaticVariants(reportGermline);
        List<ReportableGainLoss> gainsAndLosses = createCOLO829GainsLosses();
        List<ReportableGeneFusion> fusions = Lists.newArrayList();
        List<ReportableHomozygousDisruption> homozygousDisruptions = Lists.newArrayList();
        List<ReportableGeneDisruption> disruptions = createCOLO829Disruptions();
        List<ViralInsertion> viralInsertions = Lists.newArrayList();

        SampleReport sampleReport = createSkinMelanomaSampleReport(sampleId);

        String summaryWithoutGermline = "Melanoma sample showing:\n"
                + " - activating BRAF mutation that is associated with response to BRAF-inhibitors (in combination with a MEK-inhibitor)\n"
                + " - complete inactivation of CDKN2A, indicating potential benefit of CDK4/6 inhibitors\n"
                + " - complete inactivation/loss of PTEN likely resulting in an activation of the PI3K-AKT-mTOR pathway "
                + "and indicating potential benefit of mTOR/PI3K inhibitors\n"
                + " - high mutational burden (mutational load (ML) of 180, tumor mutation burden (TMB) of 13.6) that is "
                + "potentially associated with an increased response rate to checkpoint inhibitor immunotherapy";

        String summaryWithGermline = "Melanoma sample showing:\n"
                + " - activating BRAF mutation that is associated with response to BRAF-inhibitors (in combination with a MEK-inhibitor)\n"
                + " - complete inactivation of CDKN2A, indicating potential benefit of CDK4/6 inhibitors. The observed CDKN2A mutation is "
                + "also present in the germline of the patient. Referral to a genetic specialist should be considered.\n"
                + " - complete inactivation/loss of PTEN likely resulting in an activation of the PI3K-AKT-mTOR pathway "
                + "and indicating potential benefit of mTOR/PI3K inhibitors\n"
                + " - high mutational burden (mutational load (ML) of 180, tumor mutation burden (TMB) of 13.6) that is "
                + "potentially associated with an increased response rate to checkpoint inhibitor immunotherapy";

        String clinicalSummary;
        if (includeSummary && !reportGermline) {
            clinicalSummary = summaryWithoutGermline;
        } else if (includeSummary && reportGermline) {
            clinicalSummary = summaryWithGermline;
        } else {
            clinicalSummary = Strings.EMPTY;
        }

        return ImmutableAnalysedPatientReport.builder()
                .sampleReport(sampleReport)
                .qsFormNumber(qcForNumber)
                .impliedPurity(impliedTumorPurity)
                .hasReliablePurity(hasReliablePurity)
                .hasReliableQuality(true)
                .averageTumorPloidy(averageTumorPloidy)
                .clinicalSummary(clinicalSummary)
                .tumorSpecificEvidence(tumorLocationSpecificEvidence)
                .clinicalTrials(clinicalTrials)
                .offLabelEvidence(offLabelEvidence)
                .reportableVariants(reportableVariants)
                .microsatelliteIndelsPerMb(microsatelliteIndelsPerMb)
                .microsatelliteStatus(MicrosatelliteStatus.fromIndelsPerMb(microsatelliteIndelsPerMb))
                .tumorMutationalLoad(tumorMutationalLoad)
                .tumorMutationalLoadStatus(TumorMutationalStatus.fromLoad(tumorMutationalLoad))
                .tumorMutationalBurden(tumorMutationalBurden)
                .chordHrdValue(chordHrdValue)
                .chordHrdStatus(ChordStatus.fromHRD(chordHrdValue))
                .gainsAndLosses(gainsAndLosses)
                .geneFusions(fusions)
                .geneDisruptions(disruptions)
                .homozygousDisruptions(homozygousDisruptions)
                .viralInsertions(viralInsertions)
                .circosPath(CIRCOS_PATH)
                .comments(Optional.ofNullable(comments))
                .isCorrectedReport(correctionReport)
                .signaturePath(reportData.signaturePath())
                .logoRVAPath(reportData.logoRVAPath())
                .logoCompanyPath(reportData.logoCompanyPath())
                .build();
    }

    @NotNull
    public static AnalysedPatientReport buildAnalysisWithAllTablesFilledInAndReliablePurity(@NotNull String sampleId,
            @Nullable String comments) {
        return buildAnalysisWithAllTablesFilledIn(sampleId, comments, true, 1D);
    }

    @NotNull
    public static AnalysedPatientReport buildAnalysisWithAllTablesFilledIn(@NotNull String sampleId, @Nullable String comments,
            boolean hasReliablePurity, double impliedTumorPurity) {
        double averageTumorPloidy = 3.1;
        int tumorMutationalLoad = 182;
        double tumorMutationalBurden = 13.6;
        double microsatelliteIndelsPerMb = 0.1089;
        double chordHrdValue = 0.8;

        ReportData reportData = testReportData();

        List<EvidenceItem> tumorLocationSpecificEvidence = createCOLO829TumorSpecificEvidence();
        List<ClinicalTrial> clinicalTrials = createCOLO829ClinicalTrials();
        List<EvidenceItem> offLabelEvidence = createCOLO829OffLabelEvidence();
        List<ReportableVariant> reportableVariants = createAllSomaticVariants();
        List<ReportableGainLoss> gainsAndLosses = createCOLO829GainsLosses();
        List<ReportableGeneFusion> fusions = createTestFusions();
        List<ReportableGeneDisruption> disruptions = createCOLO829Disruptions();
        List<ViralInsertion> viralInsertions = createTestViralInsertions();
        List<ReportableHomozygousDisruption> homozygousDisruptions = createTestHomozygousDisruptions();

        SampleReport sampleReport = createSkinMelanomaSampleReport(sampleId);
        String clinicalSummary = Strings.EMPTY;

        return ImmutableAnalysedPatientReport.builder()
                .sampleReport(sampleReport)
                .qsFormNumber(QsFormNumber.FOR_209.display())
                .impliedPurity(impliedTumorPurity)
                .hasReliablePurity(hasReliablePurity)
                .hasReliableQuality(true)
                .averageTumorPloidy(averageTumorPloidy)
                .clinicalSummary(clinicalSummary)
                .tumorSpecificEvidence(tumorLocationSpecificEvidence)
                .clinicalTrials(clinicalTrials)
                .offLabelEvidence(offLabelEvidence)
                .reportableVariants(reportableVariants)
                .microsatelliteIndelsPerMb(microsatelliteIndelsPerMb)
                .microsatelliteStatus(MicrosatelliteStatus.fromIndelsPerMb(microsatelliteIndelsPerMb))
                .tumorMutationalLoad(tumorMutationalLoad)
                .tumorMutationalLoadStatus(TumorMutationalStatus.fromLoad(tumorMutationalLoad))
                .tumorMutationalBurden(tumorMutationalBurden)
                .chordHrdValue(chordHrdValue)
                .chordHrdStatus(ChordStatus.fromHRD(chordHrdValue))
                .gainsAndLosses(gainsAndLosses)
                .geneFusions(fusions)
                .geneDisruptions(disruptions)
                .homozygousDisruptions(homozygousDisruptions)
                .viralInsertions(viralInsertions)
                .circosPath(CIRCOS_PATH)
                .comments(Optional.ofNullable(comments))
                .isCorrectedReport(false)
                .signaturePath(reportData.signaturePath())
                .logoRVAPath(reportData.logoRVAPath())
                .logoCompanyPath(reportData.logoCompanyPath())
                .build();
    }

    @NotNull
    public static QCFailReport buildQCFailReport(@NotNull String sampleId, @NotNull QCFailReason reason) {
        SampleReport sampleReport = createSkinMelanomaSampleReport(sampleId);

        ReportData reportData = testReportData();
        return ImmutableQCFailReport.builder()
                .sampleReport(sampleReport)
                .reason(reason)
                .comments(Optional.empty())
                .isCorrectedReport(false)
                .signaturePath(reportData.signaturePath())
                .logoRVAPath(reportData.logoRVAPath())
                .logoCompanyPath(reportData.logoCompanyPath())
                .build();
    }

    @NotNull
    private static HospitalContactData createTestHospitalContactData() {
        return ImmutableHospitalContactData.builder()
                .hospitalPI("PI")
                .requesterName("Paul")
                .requesterEmail("paul@hartwig.com")
                .hospitalName("HMF Testing Center")
                .hospitalAddress("1000 AB AMSTERDAM")
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
                .refArrivalDate(LocalDate.parse("01-Jan-2020", DATE_FORMATTER))
                .tumorArrivalDate(LocalDate.parse("05-Jan-2020", DATE_FORMATTER))
                .shallowSeqPurityString(Lims.NOT_PERFORMED_STRING)
                .labProcedures("PREP013V23-QC037V20-SEQ008V25")
                .cohort("TEST")
                .projectName("TEST-001-002")
                .submissionId("SUBM")
                .hospitalContactData(createTestHospitalContactData())
                .hospitalPatientId("HOSP1")
                .hospitalPathologySampleId("PA1")
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
                .acronym("CLXH254X2101")
                .reference("EXT10453 (NL55506.078.15)")
                .build());
        trials.add(iclusionBuilder.event("BRAF p.Val600Glu")
                .scope(EvidenceScope.SPECIFIC)
                .acronym("COWBOY")
                .reference("EXT12301 (NL71732.091.19)")
                .build());
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
                .acronym("POLARIS")
                .reference("EXT11388 (NL69569.028.19)")
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
    private static List<ReportableVariant> createCOLO829SomaticVariants(boolean reportGermlineVariant) {
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
                .gDNA("7:140453136")
                .alleleReadCount(150)
                .totalReadCount(221)
                .alleleCopyNumber(4.08)
                .totalCopyNumber(6.0)
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
                .notifyClinicalGeneticist(reportGermlineVariant)
                .gDNA("9:21971153")
                .alleleReadCount(99)
                .totalReadCount(99)
                .alleleCopyNumber(1.99)
                .totalCopyNumber(1.99)
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
                .gDNA("5:1295228")
                .alleleReadCount(56)
                .totalReadCount(65)
                .alleleCopyNumber(1.74)
                .totalCopyNumber(2.0)
                .hotspot(Hotspot.HOTSPOT)
                .biallelic(true)
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
                .gDNA("2:198266779")
                .alleleReadCount(74)
                .totalReadCount(111)
                .alleleCopyNumber(2.01)
                .totalCopyNumber(3.02)
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
                .gDNA("3:189604330")
                .alleleReadCount(47)
                .totalReadCount(112)
                .alleleCopyNumber(1.67)
                .totalCopyNumber(3.98)
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
                .gDNA("3:189604330")
                .alleleReadCount(48)
                .totalReadCount(103)
                .alleleCopyNumber(2.1)
                .totalCopyNumber(4.1)
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
                .gDNA("3:81627197")
                .alleleReadCount(48)
                .totalReadCount(103)
                .alleleCopyNumber(1.3)
                .totalCopyNumber(2.5)
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
                .junctionCopyNumber(0.4)
                .build();

        ReportableGeneFusion fusion2 = ImmutableReportableGeneFusion.builder()
                .geneStart("CLCN6")
                .geneTranscriptStart("ENST00000346436")
                .geneContextStart("Intron 1")
                .geneEnd("BRAF")
                .geneTranscriptEnd("ENST00000288602")
                .geneContextEnd("Intron 8")
                .junctionCopyNumber(1D)
                .build();

        return Lists.newArrayList(fusion1, fusion2);
    }

    @NotNull
    private static List<ReportableGeneDisruption> createCOLO829Disruptions() {
        ReportableGeneDisruption disruption1 = createDisruptionBuilder().location("10q23.31")
                .gene("PTEN")
                .range("Intron 5 -> Intron 6")
                .type("DEL")
                .junctionCopyNumber(2D)
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

