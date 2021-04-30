package com.hartwig.hmftools.patientreporter;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;
import java.util.Locale;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.clinical.ImmutablePatientPrimaryTumor;
import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfig;
import com.hartwig.hmftools.common.lims.hospital.HospitalContactData;
import com.hartwig.hmftools.common.lims.hospital.ImmutableHospitalContactData;
import com.hartwig.hmftools.common.protect.ImmutableProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.copynumber.ImmutableReportableGainLoss;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.structural.linx.FusionLikelihoodType;
import com.hartwig.hmftools.common.variant.structural.linx.FusionPhasedType;
import com.hartwig.hmftools.common.variant.structural.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;
import com.hartwig.hmftools.patientreporter.algo.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.algo.GenomicAnalysis;
import com.hartwig.hmftools.patientreporter.algo.ImmutableAnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.algo.ImmutableGenomicAnalysis;
import com.hartwig.hmftools.patientreporter.cuppa.ImmutableMolecularTissueOrigin;
import com.hartwig.hmftools.patientreporter.cuppa.MolecularTissueOrigin;
import com.hartwig.hmftools.patientreporter.qcfail.ImmutableQCFailReport;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReport;
import com.hartwig.hmftools.protect.viralbreakend.ImmutableViralbreakend;
import com.hartwig.hmftools.protect.viralbreakend.Viralbreakend;
import com.hartwig.hmftools.protect.linx.ImmutableReportableGeneDisruption;
import com.hartwig.hmftools.protect.linx.ImmutableReportableHomozygousDisruption;
import com.hartwig.hmftools.protect.linx.ReportableGeneDisruption;
import com.hartwig.hmftools.protect.linx.ReportableHomozygousDisruption;
import com.hartwig.hmftools.protect.purple.ImmutableReportableVariant;
import com.hartwig.hmftools.protect.purple.ReportableVariant;
import com.hartwig.hmftools.protect.purple.ReportableVariantSource;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ExampleAnalysisTestFactory {

    private static final DateTimeFormatter DATE_FORMATTER = DateTimeFormatter.ofPattern("dd-MMM-yyyy", Locale.ENGLISH);
    private static final String CIRCOS_PATH = Resources.getResource("test_run/purple/plot/sample.circos.png").getPath();

    private ExampleAnalysisTestFactory() {
    }

    @NotNull
    public static AnalysedPatientReport buildTestReport() {
        return buildWithCOLO829Data("TestSample",
                false,
                null,
                QsFormNumber.FOR_080.display(),
                true,
                1D,
                true,
                false,
                PatientReporterTestFactory.createTestCohortConfig());
    }

    @NotNull
    public static AnalysedPatientReport buildCOLO829(@NotNull String sampleId, boolean correctionReport, @Nullable String comments,
            @NotNull LimsCohortConfig limsCohortConfig) {
        return buildWithCOLO829Data(sampleId,
                correctionReport,
                comments,
                QsFormNumber.FOR_080.display(),
                true,
                1D,
                true,
                false,
                limsCohortConfig);
    }

    @NotNull
    public static AnalysedPatientReport buildWithCOLO829Data(@NotNull String sampleId, boolean correctionReport, @Nullable String comments,
            @NotNull String qcForNumber, boolean hasReliablePurity, double impliedTumorPurity, boolean includeSummary,
            boolean reportGermline, @NotNull LimsCohortConfig limsCohortConfig) {
        double averageTumorPloidy = 3.1;
        int tumorMutationalLoad = 190;
        double tumorMutationalBurden = 13.7;
        double microsatelliteIndelsPerMb = 0.12;
        double chordHrdValue = 0D;
        ChordStatus chordStatus = ChordStatus.HR_PROFICIENT;

        ReportData reportData = PatientReporterTestFactory.loadTestReportData();

        List<ProtectEvidence> tumorSpecificEvidence = createCOLO829TumorSpecificEvidence();
        List<ProtectEvidence> clinicalTrials = createCOLO829ClinicalTrials();
        List<ProtectEvidence> offLabelEvidence = createCOLO829OffLabelEvidence();
        List<ReportableVariant> reportableVariants = createCOLO829SomaticVariants(reportGermline);
        List<ReportableGainLoss> gainsAndLosses = createCOLO829GainsLosses();
        List<LinxFusion> fusions = Lists.newArrayList();
        List<ReportableHomozygousDisruption> homozygousDisruptions = Lists.newArrayList();
        List<ReportableGeneDisruption> disruptions = createCOLO829Disruptions();
        List<Viralbreakend> viralBreakends = Lists.newArrayList();

        SampleReport sampleReport = createSkinMelanomaSampleReport(sampleId, reportGermline, limsCohortConfig);

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

        GenomicAnalysis analysis = ImmutableGenomicAnalysis.builder()
                .impliedPurity(impliedTumorPurity)
                .hasReliablePurity(hasReliablePurity)
                .hasReliableQuality(true)
                .averageTumorPloidy(averageTumorPloidy)
                .tumorSpecificEvidence(tumorSpecificEvidence)
                .clinicalTrials(clinicalTrials)
                .offLabelEvidence(offLabelEvidence)
                .reportableVariants(reportableVariants)
                .microsatelliteIndelsPerMb(microsatelliteIndelsPerMb)
                .microsatelliteStatus(MicrosatelliteStatus.fromIndelsPerMb(microsatelliteIndelsPerMb))
                .tumorMutationalLoad(tumorMutationalLoad)
                .tumorMutationalLoadStatus(TumorMutationalStatus.fromLoad(tumorMutationalLoad))
                .tumorMutationalBurden(tumorMutationalBurden)
                .chordHrdValue(chordHrdValue)
                .chordHrdStatus(chordStatus)
                .gainsAndLosses(gainsAndLosses)
                .geneFusions(fusions)
                .geneDisruptions(disruptions)
                .homozygousDisruptions(homozygousDisruptions)
                .viralBreakends(viralBreakends)
                .build();

        MolecularTissueOrigin molecularTissueOrigin =
                ImmutableMolecularTissueOrigin.builder().molecularTissueOriginResult("Skin").molecularTissueOriginPlot(CIRCOS_PATH).build();

        return ImmutableAnalysedPatientReport.builder()
                .sampleReport(sampleReport)
                .qsFormNumber(qcForNumber)
                .clinicalSummary(clinicalSummary)
                .genomicAnalysis(analysis)
                .circosPath(CIRCOS_PATH)
                .molecularTissueOrigin(molecularTissueOrigin)
                .comments(Optional.ofNullable(comments))
                .isCorrectedReport(correctionReport)
                .signaturePath(reportData.signaturePath())
                .logoRVAPath(reportData.logoRVAPath())
                .logoCompanyPath(reportData.logoCompanyPath())
                .pipelineVersion("5.19")
                .build();
    }

    @NotNull
    public static AnalysedPatientReport buildAnalysisWithAllTablesFilledInAndReliablePurity(@NotNull String sampleId,
            @Nullable String comments, @NotNull LimsCohortConfig limsCohortConfig) {
        return buildAnalysisWithAllTablesFilledIn(sampleId, comments, true, 1D, limsCohortConfig);
    }

    @NotNull
    public static AnalysedPatientReport buildAnalysisWithAllTablesFilledIn(@NotNull String sampleId, @Nullable String comments,
            boolean hasReliablePurity, double impliedTumorPurity, @NotNull LimsCohortConfig limsCohortConfig) {
        double averageTumorPloidy = 3.1;
        int tumorMutationalLoad = 182;
        double tumorMutationalBurden = 13.6;
        double microsatelliteIndelsPerMb = 0.1089;
        double chordHrdValue = 0.8;
        ChordStatus chordStatus = ChordStatus.HR_DEFICIENT;

        ReportData reportData = PatientReporterTestFactory.loadTestReportData();

        List<ProtectEvidence> tumorSpecificEvidence = createCOLO829TumorSpecificEvidence();
        List<ProtectEvidence> clinicalTrials = createCOLO829ClinicalTrials();
        List<ProtectEvidence> offLabelEvidence = createCOLO829OffLabelEvidence();
        List<ReportableVariant> reportableVariants = createAllSomaticVariants();
        List<ReportableGainLoss> gainsAndLosses = createCOLO829GainsLosses();
        List<LinxFusion> fusions = createTestFusions();
        List<ReportableGeneDisruption> disruptions = createCOLO829Disruptions();
        List<Viralbreakend> viralBreakends = createTestViralBreakends();
        List<ReportableHomozygousDisruption> homozygousDisruptions = createTestHomozygousDisruptions();

        SampleReport sampleReport = createSkinMelanomaSampleReport(sampleId, true, limsCohortConfig);
        String clinicalSummary = Strings.EMPTY;

        GenomicAnalysis analysis = ImmutableGenomicAnalysis.builder()
                .impliedPurity(impliedTumorPurity)
                .hasReliablePurity(hasReliablePurity)
                .hasReliableQuality(true)
                .averageTumorPloidy(averageTumorPloidy)
                .tumorSpecificEvidence(tumorSpecificEvidence)
                .clinicalTrials(clinicalTrials)
                .offLabelEvidence(offLabelEvidence)
                .reportableVariants(reportableVariants)
                .microsatelliteIndelsPerMb(microsatelliteIndelsPerMb)
                .microsatelliteStatus(MicrosatelliteStatus.fromIndelsPerMb(microsatelliteIndelsPerMb))
                .tumorMutationalLoad(tumorMutationalLoad)
                .tumorMutationalLoadStatus(TumorMutationalStatus.fromLoad(tumorMutationalLoad))
                .tumorMutationalBurden(tumorMutationalBurden)
                .chordHrdValue(chordHrdValue)
                .chordHrdStatus(chordStatus)
                .gainsAndLosses(gainsAndLosses)
                .geneFusions(fusions)
                .geneDisruptions(disruptions)
                .homozygousDisruptions(homozygousDisruptions)
                .viralBreakends(viralBreakends)
                .build();

        MolecularTissueOrigin molecularTissueOrigin =
                ImmutableMolecularTissueOrigin.builder().molecularTissueOriginResult("Skin").molecularTissueOriginPlot(CIRCOS_PATH).build();

        return ImmutableAnalysedPatientReport.builder()
                .sampleReport(sampleReport)
                .qsFormNumber(QsFormNumber.FOR_209.display())
                .clinicalSummary(clinicalSummary)
                .genomicAnalysis(analysis)
                .circosPath(CIRCOS_PATH)
                .molecularTissueOrigin(molecularTissueOrigin)
                .comments(Optional.ofNullable(comments))
                .isCorrectedReport(false)
                .signaturePath(reportData.signaturePath())
                .logoRVAPath(reportData.logoRVAPath())
                .logoCompanyPath(reportData.logoCompanyPath())
                .pipelineVersion("5.19")
                .build();
    }

    @NotNull
    public static QCFailReport buildQCFailReport(@NotNull String sampleId, @NotNull QCFailReason reason,
            @NotNull LimsCohortConfig limsCohortConfig) {
        SampleReport sampleReport = createSkinMelanomaSampleReport(sampleId, true, limsCohortConfig);

        ReportData reportData = PatientReporterTestFactory.loadTestReportData();
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
    private static SampleReport createSkinMelanomaSampleReport(@NotNull String sample, boolean reportGermline,
            @NotNull LimsCohortConfig cohort) {
        SampleMetadata sampleMetadata = ImmutableSampleMetadata.builder()
                .patientId("COLO829")
                .refSampleId(Strings.EMPTY)
                .refSampleBarcode("FR12123488")
                .tumorSampleId(sample)
                .tumorSampleBarcode("FR12345678")
                .build();

        return ImmutableSampleReport.builder()
                .sampleMetadata(sampleMetadata)
                .patientPrimaryTumor(ImmutablePatientPrimaryTumor.builder()
                        .patientIdentifier(sample)
                        .location("Skin")
                        .subLocation(Strings.EMPTY)
                        .type("Melanoma")
                        .subType(Strings.EMPTY)
                        .extraDetails(Strings.EMPTY)
                        .doids(Lists.newArrayList("8923"))
                        .isOverridden(false)
                        .build())
                .germlineReportingLevel(reportGermline
                        ? LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION
                        : LimsGermlineReportingLevel.NO_REPORTING)
                .reportViralInsertions(false)
                .refArrivalDate(LocalDate.parse("01-Oct-2020", DATE_FORMATTER))
                .tumorArrivalDate(LocalDate.parse("05-Oct-2020", DATE_FORMATTER))
                .shallowSeqPurityString(Lims.NOT_PERFORMED_STRING)
                .labProcedures("PREP013V23-QC037V20-SEQ008V25")
                .cohort(cohort)
                .projectName("TEST-001-002")
                .submissionId("SUBM")
                .hospitalContactData(createTestHospitalContactData())
                .hospitalPatientId("HOSP1")
                .hospitalPathologySampleId("PA1")
                .build();
    }

    @NotNull
    private static List<ProtectEvidence> createCOLO829TumorSpecificEvidence() {
        List<ProtectEvidence> evidenceItemsOnLabel = Lists.newArrayList();

        ImmutableProtectEvidence.Builder onLabelBuilder = ImmutableProtectEvidence.builder();

        evidenceItemsOnLabel.add(onLabelBuilder.genomicEvent("BRAF p.Val600Glu")
                .germline(false)
                .reported(true)
                .treatment("Cobimetinib + Vemurafenib")
                .onLabel(true)
                .level(EvidenceLevel.A)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CGI))
                .urls(Lists.newArrayList("https://www.google.com/#q=FDA"))
                .build());

        evidenceItemsOnLabel.add(onLabelBuilder.genomicEvent("BRAF p.Val600Glu")
                .germline(false)
                .reported(true)
                .treatment("Dabrafenib")
                .onLabel(true)
                .level(EvidenceLevel.A)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CGI))
                .urls(Lists.newArrayList("https://www.google.com/#q=FDA"))
                .build());

        evidenceItemsOnLabel.add(onLabelBuilder.genomicEvent("BRAF p.Val600Glu")
                .germline(false)
                .reported(true)
                .treatment("Dabrafenib + Trametinib")
                .onLabel(true)
                .level(EvidenceLevel.A)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CGI, Knowledgebase.VICC_CIVIC))
                .urls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/25399551", "https://www.google.com/#q=FDA"))
                .build());

        evidenceItemsOnLabel.add(onLabelBuilder.genomicEvent("BRAF p.Val600Glu")
                .germline(false)
                .reported(true)
                .treatment("Trametinib")
                .onLabel(true)
                .level(EvidenceLevel.A)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CGI))
                .urls(Lists.newArrayList("https://www.google.com/#q=FDA"))
                .build());

        evidenceItemsOnLabel.add(onLabelBuilder.genomicEvent("BRAF p.Val600Glu")
                .germline(false)
                .reported(true)
                .treatment("Vemurafenib")
                .onLabel(true)
                .level(EvidenceLevel.A)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CGI, Knowledgebase.VICC_CIVIC))
                .urls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/21639808",
                        "https://www.google.com/#q=FDA",
                        "https://www.google.com/#q=NCCN"))
                .build());

        evidenceItemsOnLabel.add(onLabelBuilder.genomicEvent("BRAF p.Val600Glu")
                .germline(false)
                .reported(true)
                .treatment("RO4987655")
                .onLabel(true)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CIVIC))
                .urls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/24947927"))
                .build());

        evidenceItemsOnLabel.add(onLabelBuilder.genomicEvent("PTEN partial loss")
                .germline(false)
                .reported(true)
                .treatment("Buparlisib + Carboplatin + Paclitaxel")
                .onLabel(true)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CIVIC))
                .urls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/25672916"))
                .build());

        return evidenceItemsOnLabel;
    }

    @NotNull
    private static List<ProtectEvidence> createCOLO829ClinicalTrials() {
        List<ProtectEvidence> trialsOnLabel = Lists.newArrayList();
        ImmutableProtectEvidence.Builder onLabelBuilder = ImmutableProtectEvidence.builder();

        trialsOnLabel.add(onLabelBuilder.genomicEvent("BRAF p.Val600Glu")
                .germline(false)
                .reported(true)
                .treatment("Array 818-103")
                .onLabel(true)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(Sets.newHashSet(Knowledgebase.ICLUSION))
                .urls(Lists.newArrayList("https://iclusion.org/hmf/13054"))
                .build());

        trialsOnLabel.add(onLabelBuilder.genomicEvent("BRAF p.Val600Glu")
                .germline(false)
                .reported(true)
                .treatment("CLXH254C12201")
                .onLabel(true)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(Sets.newHashSet(Knowledgebase.ICLUSION))
                .urls(Lists.newArrayList("https://iclusion.org/hmf/13660"))
                .build());

        trialsOnLabel.add(onLabelBuilder.genomicEvent("BRAF p.Val600Glu")
                .germline(false)
                .reported(true)
                .treatment("COWBOY")
                .onLabel(true)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(Sets.newHashSet(Knowledgebase.ICLUSION))
                .urls(Lists.newArrayList("https://iclusion.org/hmf/12301"))
                .build());

        trialsOnLabel.add(onLabelBuilder.genomicEvent("BRAF p.Val600Glu")
                .germline(false)
                .reported(true)
                .treatment("DRUP")
                .onLabel(true)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(Sets.newHashSet(Knowledgebase.ICLUSION))
                .urls(Lists.newArrayList("https://iclusion.org/hmf/10299"))
                .build());

        trialsOnLabel.add(onLabelBuilder.genomicEvent("BRAF p.Val600Glu")
                .germline(false)
                .reported(true)
                .treatment("EBIN (EORTC-1612-MG)")
                .onLabel(true)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(Sets.newHashSet(Knowledgebase.ICLUSION))
                .urls(Lists.newArrayList("https://iclusion.org/hmf/11284"))
                .build());

        trialsOnLabel.add(onLabelBuilder.genomicEvent("BRAF p.Val600Glu")
                .germline(false)
                .reported(true)
                .treatment("POLARIS")
                .onLabel(true)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(Sets.newHashSet(Knowledgebase.ICLUSION))
                .urls(Lists.newArrayList("https://iclusion.org/hmf/11388"))
                .build());

        trialsOnLabel.add(onLabelBuilder.genomicEvent("CDKN2A p.Ala68fs")
                .germline(false)
                .reported(true)
                .treatment("DRUP")
                .onLabel(true)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(Sets.newHashSet(Knowledgebase.ICLUSION))
                .urls(Lists.newArrayList("https://iclusion.org/hmf/10299"))
                .build());

        trialsOnLabel.add(onLabelBuilder.genomicEvent("High tumor mutation load")
                .germline(false)
                .reported(true)
                .treatment("BASKET OF BASKETS (VHIO17002)")
                .onLabel(true)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(Sets.newHashSet(Knowledgebase.ICLUSION))
                .urls(Lists.newArrayList("https://iclusion.org/hmf/11087"))
                .build());

        trialsOnLabel.add(onLabelBuilder.genomicEvent("High tumor mutation load")
                .germline(false)
                .reported(true)
                .treatment("CheckMate 848")
                .onLabel(true)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(Sets.newHashSet(Knowledgebase.ICLUSION))
                .urls(Lists.newArrayList("https://iclusion.org/hmf/10560"))
                .build());

        trialsOnLabel.add(onLabelBuilder.genomicEvent("High tumor mutation load")
                .germline(false)
                .reported(true)
                .treatment("DRUP")
                .onLabel(true)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(Sets.newHashSet(Knowledgebase.ICLUSION))
                .urls(Lists.newArrayList("https://iclusion.org/hmf/10299"))
                .build());

        trialsOnLabel.add(onLabelBuilder.genomicEvent("PTEN partial loss")
                .germline(false)
                .reported(true)
                .treatment("DRUP")
                .onLabel(true)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(Sets.newHashSet(Knowledgebase.ICLUSION))
                .urls(Lists.newArrayList("https://iclusion.org/hmf/10299"))
                .build());

        return trialsOnLabel;
    }

    @NotNull
    private static List<ProtectEvidence> createCOLO829OffLabelEvidence() {
        List<ProtectEvidence> evidenceItemsOffLabel = Lists.newArrayList();

        ImmutableProtectEvidence.Builder onLabelBuilder = ImmutableProtectEvidence.builder();

        evidenceItemsOffLabel.add(onLabelBuilder.genomicEvent("BRAF p.Val600Glu")
                .germline(false)
                .reported(true)
                .treatment("Bevacizumab")
                .onLabel(false)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CIVIC))
                .urls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/19571295", "http://www.ncbi.nlm.nih.gov/pubmed/19603024"))
                .build());

        evidenceItemsOffLabel.add(onLabelBuilder.genomicEvent("BRAF p.Val600Glu")
                .germline(false)
                .reported(true)
                .treatment("CI-1040")
                .onLabel(false)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CIVIC))
                .urls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/18682506", "http://www.ncbi.nlm.nih.gov/pubmed/21882184"))
                .build());

        evidenceItemsOffLabel.add(onLabelBuilder.genomicEvent("BRAF p.Val600Glu")
                .germline(false)
                .reported(true)
                .treatment("Cetuximab")
                .onLabel(false)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CIVIC))
                .urls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/25673558"))
                .build());

        evidenceItemsOffLabel.add(onLabelBuilder.genomicEvent("BRAF p.Val600Glu")
                .germline(false)
                .reported(true)
                .treatment("Cetuximab")
                .onLabel(false)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CGI, Knowledgebase.VICC_CIVIC))
                .urls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/19001320",
                        "http://www.ncbi.nlm.nih.gov/pubmed/19884556",
                        "http://www.ncbi.nlm.nih.gov/pubmed/20619739",
                        "http://www.ncbi.nlm.nih.gov/pubmed/21163703",
                        "http://www.ncbi.nlm.nih.gov/pubmed/23325582",
                        "http://www.ncbi.nlm.nih.gov/pubmed/25666295",
                        "http://www.ncbi.nlm.nih.gov/pubmed/25989278"))
                .build());

        evidenceItemsOffLabel.add(onLabelBuilder.genomicEvent("BRAF p.Val600Glu")
                .germline(false)
                .reported(true)
                .treatment("Cetuximab + Irinotecan + Vemurafenib")
                .onLabel(false)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CIVIC))
                .urls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/27729313"))
                .build());

        evidenceItemsOffLabel.add(onLabelBuilder.genomicEvent("BRAF p.Val600Glu")
                .germline(false)
                .reported(true)
                .treatment("Cetuximab + Vemurafenib")
                .onLabel(false)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CIVIC))
                .urls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/26287849"))
                .build());

        evidenceItemsOffLabel.add(onLabelBuilder.genomicEvent("BRAF p.Val600Glu")
                .germline(false)
                .reported(true)
                .treatment("Fluorouracil")
                .onLabel(false)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CIVIC))
                .urls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/19603024"))
                .build());

        evidenceItemsOffLabel.add(onLabelBuilder.genomicEvent("BRAF p.Val600Glu")
                .germline(false)
                .reported(true)
                .treatment("Irinotecan")
                .onLabel(false)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CIVIC))
                .urls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/19603024"))
                .build());

        evidenceItemsOffLabel.add(onLabelBuilder.genomicEvent("BRAF p.Val600Glu")
                .germline(false)
                .reported(true)
                .treatment("Oxaliplatin")
                .onLabel(false)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CIVIC))
                .urls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/19603024"))
                .build());

        evidenceItemsOffLabel.add(onLabelBuilder.genomicEvent("BRAF p.Val600Glu")
                .germline(false)
                .reported(true)
                .treatment("Panitumumab")
                .onLabel(false)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CIVIC))
                .urls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/25673558"))
                .build());

        evidenceItemsOffLabel.add(onLabelBuilder.genomicEvent("BRAF p.Val600Glu")
                .germline(false)
                .reported(true)
                .treatment("Panitumumab")
                .onLabel(false)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CGI, Knowledgebase.VICC_CIVIC))
                .urls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/19001320",
                        "http://www.ncbi.nlm.nih.gov/pubmed/20619739",
                        "http://www.ncbi.nlm.nih.gov/pubmed/21163703",
                        "http://www.ncbi.nlm.nih.gov/pubmed/23325582",
                        "http://www.ncbi.nlm.nih.gov/pubmed/25989278"))
                .build());

        evidenceItemsOffLabel.add(onLabelBuilder.genomicEvent("BRAF p.Val600Glu")
                .germline(false)
                .reported(true)
                .treatment("Panitumumab + Vemurafenib")
                .onLabel(false)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CIVIC))
                .urls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/25589621"))
                .build());

        evidenceItemsOffLabel.add(onLabelBuilder.genomicEvent("BRAF p.Val600Glu")
                .germline(false)
                .reported(true)
                .treatment("Selumetinib")
                .onLabel(false)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CIVIC))
                .urls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/22492957"))
                .build());

        evidenceItemsOffLabel.add(onLabelBuilder.genomicEvent("BRAF p.Val600Glu")
                .germline(false)
                .reported(true)
                .treatment("Sorafenib")
                .onLabel(false)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CIVIC))
                .urls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/18682506", "http://www.ncbi.nlm.nih.gov/pubmed/21882184"))
                .build());

        evidenceItemsOffLabel.add(onLabelBuilder.genomicEvent("BRAF p.Val600Glu")
                .germline(false)
                .reported(true)
                .treatment("Vemurafenib")
                .onLabel(false)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CIVIC))
                .urls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/26287849"))
                .build());

        evidenceItemsOffLabel.add(onLabelBuilder.genomicEvent("PTEN partial loss")
                .germline(false)
                .reported(true)
                .treatment("Anti-EGFR monoclonal antibody")
                .onLabel(false)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CGI))
                .urls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/19398573", "http://www.ncbi.nlm.nih.gov/pubmed/21163703"))
                .build());

        evidenceItemsOffLabel.add(onLabelBuilder.genomicEvent("PTEN partial loss")
                .germline(false)
                .reported(true)
                .treatment("Cetuximab")
                .onLabel(false)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CIVIC))
                .urls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/21163703"))
                .build());

        evidenceItemsOffLabel.add(onLabelBuilder.genomicEvent("PTEN partial loss")
                .germline(false)
                .reported(true)
                .treatment("Everolimus")
                .onLabel(false)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CIVIC))
                .urls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/23989949"))
                .build());

        evidenceItemsOffLabel.add(onLabelBuilder.genomicEvent("PTEN partial loss")
                .germline(false)
                .reported(true)
                .treatment("Lapatinib + Trastuzumab")
                .onLabel(false)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CIVIC))
                .urls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/25300346"))
                .build());

        evidenceItemsOffLabel.add(onLabelBuilder.genomicEvent("PTEN partial loss")
                .germline(false)
                .reported(true)
                .treatment("Ridaforolimus")
                .onLabel(false)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CIVIC))
                .urls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/24166148"))
                .build());

        evidenceItemsOffLabel.add(onLabelBuilder.genomicEvent("PTEN partial loss")
                .germline(false)
                .reported(true)
                .treatment("Temsirolimus")
                .onLabel(false)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CIVIC))
                .urls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/24166148"))
                .build());

        evidenceItemsOffLabel.add(onLabelBuilder.genomicEvent("PTEN partial loss")
                .germline(false)
                .reported(true)
                .treatment("Trastuzumab")
                .onLabel(false)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CIVIC))
                .urls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/20813970", "http://www.ncbi.nlm.nih.gov/pubmed/24387334"))
                .build());

        return evidenceItemsOffLabel;
    }

    @NotNull
    private static List<ReportableVariant> createCOLO829SomaticVariants(boolean forceCDKN2AVariantToBeGermline) {
        ReportableVariant variant1 = ImmutableReportableVariant.builder()
                .source(ReportableVariantSource.SOMATIC)
                .gene("BRAF")
                .genotypeStatus(GenotypeStatus.HET)
                .chromosome("7")
                .position(140453136)
                .ref("T")
                .alt("A")
                .type(VariantType.SNP)
                .canonicalTranscript("ENST00000288602")
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .canonicalHgvsCodingImpact("c.1799T>A")
                .canonicalHgvsProteinImpact("p.Val600Glu")
                .alleleReadCount(150)
                .totalReadCount(221)
                .alleleCopyNumber(4.08)
                .totalCopyNumber(6.0)
                .hotspot(Hotspot.HOTSPOT)
                .driverLikelihood(1D)
                .clonalLikelihood(1D)
                .biallelic(false)
                .build();

        ReportableVariant variant2 = ImmutableReportableVariant.builder()
                .source(forceCDKN2AVariantToBeGermline ? ReportableVariantSource.GERMLINE : ReportableVariantSource.SOMATIC)
                .gene("CDKN2A")
                .genotypeStatus(GenotypeStatus.HET)
                .chromosome("9")
                .position(21971153)
                .ref("CCG")
                .alt("C")
                .type(VariantType.INDEL)
                .canonicalTranscript("ENST00000498124")
                .canonicalCodingEffect(CodingEffect.NONSENSE_OR_FRAMESHIFT)
                .canonicalHgvsCodingImpact("c.203_204delCG")
                .canonicalHgvsProteinImpact("p.Ala68fs")
                .alleleReadCount(99)
                .totalReadCount(99)
                .alleleCopyNumber(1.99)
                .totalCopyNumber(1.99)
                .hotspot(Hotspot.NEAR_HOTSPOT)
                .clonalLikelihood(1D)
                .driverLikelihood(0.9)
                .biallelic(true)
                .build();

        ReportableVariant variant3 = ImmutableReportableVariant.builder()
                .source(ReportableVariantSource.SOMATIC)
                .gene("TERT")
                .genotypeStatus(GenotypeStatus.HET)
                .chromosome("5")
                .position(1295228)
                .ref("GG")
                .alt("AA")
                .type(VariantType.MNP)
                .canonicalTranscript("ENST00000310581")
                .canonicalCodingEffect(CodingEffect.NONE)
                .canonicalHgvsCodingImpact("c.-125_-124delCCinsTT")
                .canonicalHgvsProteinImpact(Strings.EMPTY)
                .alleleReadCount(56)
                .totalReadCount(65)
                .alleleCopyNumber(1.74)
                .totalCopyNumber(2.0)
                .hotspot(Hotspot.HOTSPOT)
                .clonalLikelihood(1D)
                .driverLikelihood(0.85)
                .biallelic(true)
                .build();

        ReportableVariant variant4 = ImmutableReportableVariant.builder()
                .source(ReportableVariantSource.SOMATIC)
                .gene("SF3B1")
                .genotypeStatus(GenotypeStatus.HET)
                .chromosome("2")
                .position(198266779)
                .ref("C")
                .alt("T")
                .type(VariantType.SNP)
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .canonicalTranscript("ENST00000335508")
                .canonicalHgvsCodingImpact("c.2153C>T")
                .canonicalHgvsProteinImpact("p.Pro718Leu")
                .alleleReadCount(74)
                .totalReadCount(111)
                .alleleCopyNumber(2.01)
                .totalCopyNumber(3.02)
                .hotspot(Hotspot.NON_HOTSPOT)
                .clonalLikelihood(1D)
                .driverLikelihood(0.15)
                .biallelic(false)
                .build();

        ReportableVariant variant5 = ImmutableReportableVariant.builder()
                .source(ReportableVariantSource.SOMATIC)
                .gene("TP63")
                .genotypeStatus(GenotypeStatus.HET)
                .chromosome("3")
                .position(189604330)
                .ref("G")
                .alt("T")
                .type(VariantType.SNP)
                .canonicalTranscript("ENST00000264731")
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .canonicalHgvsCodingImpact("c.1497G>T")
                .canonicalHgvsProteinImpact("p.Met499Ile")
                .alleleReadCount(47)
                .totalReadCount(112)
                .alleleCopyNumber(1.67)
                .totalCopyNumber(3.98)
                .hotspot(Hotspot.NON_HOTSPOT)
                .clonalLikelihood(1D)
                .driverLikelihood(0.1)
                .biallelic(false)
                .build();

        return Lists.newArrayList(variant1, variant2, variant3, variant4, variant5);
    }

    @NotNull
    private static List<ReportableVariant> createAllSomaticVariants() {
        ReportableVariant variant1 = ImmutableReportableVariant.builder()
                .source(ReportableVariantSource.SOMATIC)
                .gene("TP63")
                .genotypeStatus(GenotypeStatus.HET)
                .chromosome("3")
                .position(189604330)
                .ref("G")
                .alt("T")
                .type(VariantType.SNP)
                .canonicalTranscript("123")
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .canonicalHgvsCodingImpact("c.1497G>T")
                .canonicalHgvsProteinImpact("p.Met499Ile")
                .alleleReadCount(48)
                .totalReadCount(103)
                .alleleCopyNumber(2.1)
                .totalCopyNumber(4.1)
                .hotspot(Hotspot.NON_HOTSPOT)
                .clonalLikelihood(0.47)
                .driverLikelihood(0.1)
                .biallelic(false)
                .build();

        ReportableVariant variant2 = ImmutableReportableVariant.builder()
                .source(ReportableVariantSource.SOMATIC)
                .gene("KIT")
                .genotypeStatus(GenotypeStatus.HET)
                .chromosome("3")
                .position(81627197)
                .ref("G")
                .alt("T")
                .type(VariantType.SNP)
                .canonicalTranscript("123")
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .canonicalHgvsCodingImpact("c.1497G>T")
                .canonicalHgvsProteinImpact("p.Met499Ile")
                .alleleReadCount(48)
                .totalReadCount(103)
                .alleleCopyNumber(1.3)
                .totalCopyNumber(2.5)
                .hotspot(Hotspot.NON_HOTSPOT)
                .clonalLikelihood(0.68)
                .driverLikelihood(0.1)
                .biallelic(true)
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
    private static List<LinxFusion> createTestFusions() {
        LinxFusion fusion1 = ImmutableLinxFusion.builder()
                .fivePrimeBreakendId(1)
                .threePrimeBreakendId(2)
                .name(Strings.EMPTY)
                .reported(true)
                .reportedType(Strings.EMPTY)
                .phased(FusionPhasedType.INFRAME)
                .likelihood(FusionLikelihoodType.HIGH)
                .chainLength(1)
                .chainLinks(1)
                .chainTerminated(true)
                .domainsKept(Strings.EMPTY)
                .domainsLost(Strings.EMPTY)
                .skippedExonsUp(2)
                .skippedExonsDown(4)
                .fusedExonUp(6)
                .fusedExonDown(7)
                .geneStart("TMPRSS2")
                .geneContextStart("Intron 5")
                .geneTranscriptStart("ENST00000398585")
                .geneEnd("PNPLA7")
                .geneContextEnd("Intron 3")
                .geneTranscriptEnd("ENST00000406427")
                .junctionCopyNumber(0.4)
                .build();

        LinxFusion fusion2 = ImmutableLinxFusion.builder()
                .fivePrimeBreakendId(1)
                .threePrimeBreakendId(2)
                .name(Strings.EMPTY)
                .reported(true)
                .reportedType(Strings.EMPTY)
                .phased(FusionPhasedType.SKIPPED_EXONS)
                .likelihood(FusionLikelihoodType.LOW)
                .chainLength(1)
                .chainLinks(1)
                .chainTerminated(true)
                .domainsKept(Strings.EMPTY)
                .domainsLost(Strings.EMPTY)
                .skippedExonsUp(2)
                .skippedExonsDown(4)
                .fusedExonUp(6)
                .fusedExonDown(7)
                .geneStart("CLCN6")
                .geneContextStart("Intron 1")
                .geneTranscriptStart("ENST00000346436")
                .geneEnd("BRAF")
                .geneContextEnd("Intron 8")
                .geneTranscriptEnd("ENST00000288602")
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
                .firstAffectedExon(5)
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
    private static List<Viralbreakend> createTestViralBreakends() {
        List<Viralbreakend> viralbreakends = Lists.newArrayList(ImmutableViralbreakend.builder()
                .taxidGenus(Strings.EMPTY)
                .nameGenus(Strings.EMPTY)
                .readsGenusTree(Strings.EMPTY)
                .taxidSpecies(Strings.EMPTY)
                .nameSpecies(Strings.EMPTY)
                .readsSpeciesTree(Strings.EMPTY)
                .taxidAssigned(Strings.EMPTY)
                .nameAssigned("Human papillomavirus type 16")
                .readsAssignedTree(Strings.EMPTY)
                .readsAssignedDirect(Strings.EMPTY)
                .Reference(Strings.EMPTY)
                .referenceTaxid(Strings.EMPTY)
                .referenceKmerCount(Strings.EMPTY)
                .alternateKmerCountRname(Strings.EMPTY)
                .startpos(Strings.EMPTY)
                .endpos(Strings.EMPTY)
                .numreads(Strings.EMPTY)
                .covbases(Strings.EMPTY)
                .coverage(Strings.EMPTY)
                .meandepth(Strings.EMPTY)
                .meanbaseq(Strings.EMPTY)
                .meanmapq(Strings.EMPTY)
                .integrations("2")
                .QCStatus(Strings.EMPTY)
                .build());
        return Lists.newArrayList(viralbreakends);
    }

    @NotNull
    private static ImmutableReportableGeneDisruption.Builder createDisruptionBuilder() {
        return ImmutableReportableGeneDisruption.builder().firstAffectedExon(1);
    }
}

