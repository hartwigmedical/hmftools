package com.hartwig.hmftools.patientreporter;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Optional;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.clinical.ImmutablePatientPrimaryTumor;
import com.hartwig.hmftools.common.cuppa.ImmutableMolecularTissueOrigin;
import com.hartwig.hmftools.common.cuppa.MolecularTissueOrigin;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfig;
import com.hartwig.hmftools.common.lims.hospital.HospitalContactData;
import com.hartwig.hmftools.common.lims.hospital.ImmutableHospitalContactData;
import com.hartwig.hmftools.common.linx.ImmutableReportableGeneDisruption;
import com.hartwig.hmftools.common.linx.ImmutableReportableHomozygousDisruption;
import com.hartwig.hmftools.common.linx.ReportableGeneDisruption;
import com.hartwig.hmftools.common.linx.ReportableHomozygousDisruption;
import com.hartwig.hmftools.common.peach.ImmutablePeachGenotype;
import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.common.protect.ImmutableProtectEvidence;
import com.hartwig.hmftools.common.protect.ImmutableProtectSource;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidenceType;
import com.hartwig.hmftools.common.protect.ProtectSource;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.purple.cnchromosome.CnPerChromosomeArmData;
import com.hartwig.hmftools.common.purple.cnchromosome.ImmutableCnPerChromosomeArmData;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.copynumber.ImmutableReportableGainLoss;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.purple.segment.ChromosomeArm;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.common.sv.linx.FusionLikelihoodType;
import com.hartwig.hmftools.common.sv.linx.FusionPhasedType;
import com.hartwig.hmftools.common.sv.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.common.utils.DataUtil;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.ImmutableReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariantSource;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.VirusTestFactory;
import com.hartwig.hmftools.patientreporter.algo.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.algo.GenomicAnalysis;
import com.hartwig.hmftools.patientreporter.algo.ImmutableAnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.algo.ImmutableGenomicAnalysis;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ExampleAnalysisTestFactory {

    private static final DateTimeFormatter DATE_FORMATTER = DateTimeFormatter.ofPattern("dd-MMM-yyyy", Locale.ENGLISH);
    private static final PatientReporterConfig REPORTER_CONFIG = PatientReporterTestFactory.createTestReporterConfig();

    private static final String UDI_DI = "(01) 8720299486010(8012)v5.25";

    private ExampleAnalysisTestFactory() {
    }

    @NotNull
    public static AnalysedPatientReport createTestReport() {
        return createWithCOLO829Data(new ExampleAnalysisConfig.Builder().build(), PurpleQCStatus.PASS);
    }

    @NotNull
    public static AnalysedPatientReport createWithCOLO829Data(@NotNull ExampleAnalysisConfig config,
            @NotNull PurpleQCStatus purpleQCStatus) {
        String pipelineVersion = "5.25";
        double averageTumorPloidy = 3.1;
        int tumorMutationalLoad = 189;
        double tumorMutationalBurden = 13.7;
        double microsatelliteIndelsPerMb = 0.12;
        double chordHrdValue = 0D;
        ChordStatus chordStatus = ChordStatus.HR_PROFICIENT;
        String reportDate = DataUtil.formatDate(LocalDate.now());

        ReportData reportData = PatientReporterTestFactory.loadTestReportData();

        List<ProtectEvidence> tumorSpecificEvidence = createCOLO829TumorSpecificEvidence();
        List<ProtectEvidence> clinicalTrials = createCOLO829ClinicalTrials();
        List<ProtectEvidence> offLabelEvidence = createCOLO829OffLabelEvidence();
        List<ReportableVariant> reportableVariants = createCOLO829SomaticVariants(config.reportGermline());
        Map<ReportableVariant, Boolean> notifyGermlineStatusPerVariant = notifyAllGermlineVariants(reportableVariants);
        List<ReportableGainLoss> gainsAndLosses = createCOLO829GainsLosses();
        List<LinxFusion> fusions = Lists.newArrayList();
        List<ReportableHomozygousDisruption> homozygousDisruptions = Lists.newArrayList();
        List<ReportableGeneDisruption> disruptions = createCOLO829Disruptions();
        List<AnnotatedVirus> viruses = Lists.newArrayList();
        List<PeachGenotype> peachGenotypes = createTestPeachGenotypes();

        SampleReport sampleReport = createSkinMelanomaSampleReport(config.sampleId(), config.reportGermline(), config.limsCohortConfig());

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
        if (config.includeSummary() && !config.reportGermline()) {
            clinicalSummary = config.reportGermline() ? summaryWithGermline : summaryWithoutGermline;
        } else {
            clinicalSummary = Strings.EMPTY;
        }

        GenomicAnalysis analysis = ImmutableGenomicAnalysis.builder()
                .purpleQCStatus(Sets.newHashSet(purpleQCStatus))
                .impliedPurity(config.impliedTumorPurity())
                .hasReliablePurity(config.hasReliablePurity())
                .hasReliableQuality(true)
                .averageTumorPloidy(averageTumorPloidy)
                .tumorSpecificEvidence(tumorSpecificEvidence)
                .clinicalTrials(clinicalTrials)
                .offLabelEvidence(offLabelEvidence)
                .reportableVariants(reportableVariants)
                .notifyGermlineStatusPerVariant(notifyGermlineStatusPerVariant)
                .microsatelliteIndelsPerMb(microsatelliteIndelsPerMb)
                .microsatelliteStatus(MicrosatelliteStatus.fromIndelsPerMb(microsatelliteIndelsPerMb))
                .tumorMutationalLoad(tumorMutationalLoad)
                .tumorMutationalLoadStatus(TumorMutationalStatus.fromLoad(tumorMutationalLoad))
                .tumorMutationalBurden(tumorMutationalBurden)
                .chordHrdValue(chordHrdValue)
                .chordHrdStatus(chordStatus)
                .gainsAndLosses(gainsAndLosses)
                .cnPerChromosome(extractCnPerChromosome())
                .geneFusions(fusions)
                .geneDisruptions(disruptions)
                .homozygousDisruptions(homozygousDisruptions)
                .reportableViruses(viruses)
                .build();

        MolecularTissueOrigin molecularTissueOrigin = ImmutableMolecularTissueOrigin.builder()
                .conclusion("Melanoma (likelihood=99.6%)")
                .plotPath(REPORTER_CONFIG.molecularTissueOriginPlot())
                .build();

        return ImmutableAnalysedPatientReport.builder()
                .sampleReport(sampleReport)
                .qsFormNumber(config.qcForNumber().display())
                .clinicalSummary(clinicalSummary)
                .genomicAnalysis(analysis)
                .circosPath(REPORTER_CONFIG.purpleCircosPlot())
                .molecularTissueOrigin(molecularTissueOrigin)
                .comments(Optional.ofNullable(config.comments()))
                .isCorrectedReport(config.isCorrectionReport())
                .isCorrectedReportExtern(config.isCorrectionReportExtern())
                .signaturePath(reportData.signaturePath())
                .udiDi(UDI_DI)
                .logoRVAPath(reportData.logoRVAPath())
                .logoCompanyPath(reportData.logoCompanyPath())
                .pipelineVersion(pipelineVersion)
                .peachGenotypes(peachGenotypes)
                .reportDate(reportDate)
                .isWGSreport(true)
                .build();
    }

    @NotNull
    public static AnalysedPatientReport createAnalysisWithAllTablesFilledIn(@NotNull ExampleAnalysisConfig config,
            @NotNull PurpleQCStatus purpleQCStatus) {
        AnalysedPatientReport coloReport = createWithCOLO829Data(config, purpleQCStatus);

        List<LinxFusion> fusions = createTestFusions();
        List<AnnotatedVirus> viruses = createTestAnnotatedViruses();
        List<ReportableHomozygousDisruption> homozygousDisruptions = createTestHomozygousDisruptions();

        GenomicAnalysis analysis = ImmutableGenomicAnalysis.builder()
                .from(coloReport.genomicAnalysis())
                .geneFusions(fusions)
                .homozygousDisruptions(homozygousDisruptions)
                .reportableViruses(viruses)
                .build();

        return ImmutableAnalysedPatientReport.builder().from(coloReport).genomicAnalysis(analysis).build();
    }

    @NotNull
    private static Map<ReportableVariant, Boolean> notifyAllGermlineVariants(@NotNull List<ReportableVariant> reportableVariants) {
        Map<ReportableVariant, Boolean> notifyGermlineStatusPerVariant = Maps.newHashMap();
        for (ReportableVariant variant : reportableVariants) {
            notifyGermlineStatusPerVariant.put(variant, variant.source() == ReportableVariantSource.GERMLINE);
        }
        return notifyGermlineStatusPerVariant;
    }

    @NotNull
    public static CnPerChromosomeArmData buildCnPerChromosomeArmData(@NotNull HumanChromosome chromosome,
            @NotNull ChromosomeArm chromosomeArm, double copyNumber) {
        return ImmutableCnPerChromosomeArmData.builder().chromosome(chromosome).chromosomeArm(chromosomeArm).copyNumber(copyNumber).build();
    }

    @NotNull
    public static List<CnPerChromosomeArmData> extractCnPerChromosome() {
        List<CnPerChromosomeArmData> cnPerChromosomeArm = Lists.newArrayList();
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("1"), ChromosomeArm.P_ARM, 2.5764959002046512));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("1"), ChromosomeArm.Q_ARM, 3.922154509665307));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("2"), ChromosomeArm.P_ARM, 3.0171634513146657));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("2"), ChromosomeArm.Q_ARM, 3.0219000202305364));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("3"), ChromosomeArm.P_ARM, 3.5912655243657037));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("3"), ChromosomeArm.Q_ARM, 4.000405698398538));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("4"), ChromosomeArm.P_ARM, 2.0229999604574793));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("4"), ChromosomeArm.Q_ARM, 3.8454729078639636));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("5"), ChromosomeArm.P_ARM, 2.002090592970043));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("5"), ChromosomeArm.Q_ARM, 2.011425950136734));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("6"), ChromosomeArm.P_ARM, 3.845563676541185));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("6"), ChromosomeArm.Q_ARM, 2.9144758693840416));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("7"), ChromosomeArm.P_ARM, 4.024705530627151));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("7"), ChromosomeArm.Q_ARM, 4.169394819739314));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("8"), ChromosomeArm.P_ARM, 3.33329992648033));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("8"), ChromosomeArm.Q_ARM, 3.344172929126994));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("9"), ChromosomeArm.P_ARM, 2.7299876766236433));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("9"), ChromosomeArm.Q_ARM, 3.706061264689252));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("10"), ChromosomeArm.P_ARM, 2.502865993794371));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("10"), ChromosomeArm.Q_ARM, 2.0093221707856945));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("11"), ChromosomeArm.P_ARM, 3.1662562322138417));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("11"), ChromosomeArm.Q_ARM, 2.911332199188708));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("12"), ChromosomeArm.P_ARM, 3.0119999171541836));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("12"), ChromosomeArm.Q_ARM, 3.0024718089108817));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("13"), ChromosomeArm.P_ARM, 3.157299819582857));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("13"), ChromosomeArm.Q_ARM, 3.1479621008464864));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("14"), ChromosomeArm.P_ARM, 3.03029982684));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("14"), ChromosomeArm.Q_ARM, 3.0134803904104572));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("15"), ChromosomeArm.P_ARM, 3.7027997998486484));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("15"), ChromosomeArm.Q_ARM, 2.5464224756588587));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("16"), ChromosomeArm.P_ARM, 3.187989400854891));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("16"), ChromosomeArm.Q_ARM, 1.9895662504676845));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("17"), ChromosomeArm.P_ARM, 2.988699874228875));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("17"), ChromosomeArm.Q_ARM, 3.04380005299814));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("18"), ChromosomeArm.P_ARM, 2.370029828320411));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("18"), ChromosomeArm.Q_ARM, 2.850749440994104));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("19"), ChromosomeArm.P_ARM, 2.885974468288675));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("19"), ChromosomeArm.Q_ARM, 2.9279000888664264));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("20"), ChromosomeArm.P_ARM, 4.016485962287397));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("20"), ChromosomeArm.Q_ARM, 4.00558480110238));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("21"), ChromosomeArm.P_ARM, 2.9929997659548166));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("21"), ChromosomeArm.Q_ARM, 3.0001645829865997));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("22"), ChromosomeArm.P_ARM, 3.9840997252344827));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("22"), ChromosomeArm.Q_ARM, 3.9767647179863497));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("X"), ChromosomeArm.P_ARM, 1.9504007026407164));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("X"), ChromosomeArm.Q_ARM, 1.9559000205584387));
        return cnPerChromosomeArm;
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
                .biopsyLocation("Skin")
                .germlineReportingLevel(reportGermline
                        ? LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION
                        : LimsGermlineReportingLevel.NO_REPORTING)
                .reportViralPresence(cohort.reportViral())
                .reportPharmogenetics(cohort.reportPeach())
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
    private static Set<ProtectSource> createTestProtectSource(@NotNull Set<Knowledgebase> sources,
            @NotNull ProtectEvidenceType protectEvidenceType, @Nullable Integer rank) {
        Set<ProtectSource> protectSources = Sets.newHashSet();

        for (Knowledgebase knowledgebase : sources) {
            protectSources.add(ImmutableProtectSource.builder()
                    .sources(knowledgebase)
                    .sourceEvent(Strings.EMPTY)
                    .sourceUrls(Sets.newHashSet())
                    .evidenceType(protectEvidenceType)
                    .rangeRank(rank)
                    .build());
        }
        return protectSources;
    }

    @NotNull
    private static List<ProtectEvidence> createCOLO829TumorSpecificEvidence() {
        List<ProtectEvidence> evidenceItemsOnLabel = Lists.newArrayList();

        ImmutableProtectEvidence.Builder onLabelBuilder = ImmutableProtectEvidence.builder().onLabel(true);

        evidenceItemsOnLabel.add(onLabelBuilder.gene("BRAF")
                .transcript("123")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("Cobimetinib + Vemurafenib")
                .level(EvidenceLevel.A)
                .direction(EvidenceDirection.RESPONSIVE)
                .evidenceUrls(Lists.newArrayList("https://www.google.com/#q=FDA"))
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.VICC_CIVIC, Knowledgebase.VICC_CGI),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null))
                .build());

        evidenceItemsOnLabel.add(onLabelBuilder.gene("BRAF")
                .transcript("123")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("Dabrafenib")
                .level(EvidenceLevel.A)
                .direction(EvidenceDirection.RESPONSIVE)
                .evidenceUrls(Lists.newArrayList("https://www.google.com/#q=FDA"))
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.VICC_CIVIC, Knowledgebase.VICC_CGI),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null))
                .build());

        evidenceItemsOnLabel.add(onLabelBuilder.gene("BRAF")
                .transcript("123")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("Dabrafenib + Trametinib")
                .level(EvidenceLevel.A)
                .direction(EvidenceDirection.RESPONSIVE)
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.VICC_CIVIC, Knowledgebase.VICC_CGI),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null))
                .evidenceUrls(Lists.newArrayList("https://www.google.com/#q=FDA"))
                .build());

        evidenceItemsOnLabel.add(onLabelBuilder.gene("BRAF")
                .transcript("123")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("Trametinib")
                .level(EvidenceLevel.A)
                .direction(EvidenceDirection.RESPONSIVE)
                .evidenceUrls(Lists.newArrayList("https://www.google.com/#q=FDA"))
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.VICC_CGI), ProtectEvidenceType.EXON_MUTATION, 1))
                .build());

        evidenceItemsOnLabel.add(onLabelBuilder.gene("BRAF")
                .transcript("123")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("Vemurafenib")
                .level(EvidenceLevel.A)
                .direction(EvidenceDirection.RESPONSIVE)
                .evidenceUrls(Lists.newArrayList("https://www.google.com/#q=NCCN",
                        "http://www.ncbi.nlm.nih.gov/pubmed/21639808",
                        "https://www.google.com/#q=FDA"))
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.VICC_CGI), ProtectEvidenceType.CODON_MUTATION, 600))
                .build());

        evidenceItemsOnLabel.add(onLabelBuilder.gene("BRAF")
                .transcript("123")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("RO4987655")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .evidenceUrls(Lists.newArrayList("https://www.google.com/#q=FDA"))
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.VICC_CIVIC),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null))
                .build());

        evidenceItemsOnLabel.add(onLabelBuilder.gene("PTEN")
                .transcript("123")
                .isCanonical(true)
                .event("partial loss")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("Buparlisib + Carboplatin + Paclitaxel")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .evidenceUrls(Lists.newArrayList("https://www.google.com/#q=FDA"))
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.VICC_CGI), ProtectEvidenceType.DELETION, null))
                .build());

        return evidenceItemsOnLabel;
    }

    @NotNull
    private static List<ProtectEvidence> createCOLO829ClinicalTrials() {
        List<ProtectEvidence> trialsOnLabel = Lists.newArrayList();

        ImmutableProtectEvidence.Builder onLabelBuilder = ImmutableProtectEvidence.builder().onLabel(true);

        trialsOnLabel.add(onLabelBuilder.gene("BRAF")
                .transcript("123")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("Array 818-103")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .evidenceUrls(Lists.newArrayList("https://www.google.com/#q=FDA"))
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.ICLUSION),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null))
                .build());

        trialsOnLabel.add(onLabelBuilder.gene("BRAF")
                .transcript("123")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("CLXH254C12201")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .evidenceUrls(Lists.newArrayList("https://www.google.com/#q=FDA"))
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.ICLUSION),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null))
                .build());

        trialsOnLabel.add(onLabelBuilder.gene("BRAF")
                .transcript("123")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("COWBOY")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .evidenceUrls(Lists.newArrayList("https://www.google.com/#q=FDA"))
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.ICLUSION),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null))
                .build());

        trialsOnLabel.add(onLabelBuilder.gene("BRAF")
                .transcript("123")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("COWBOY")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .evidenceUrls(Lists.newArrayList("https://www.google.com/#q=FDA"))
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.ICLUSION),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null))
                .build());

        trialsOnLabel.add(onLabelBuilder.gene("BRAF")
                .transcript("123")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("DRUP")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .evidenceUrls(Lists.newArrayList("https://www.google.com/#q=FDA"))
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.ICLUSION),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null))
                .build());

        trialsOnLabel.add(onLabelBuilder.gene("BRAF")
                .transcript("123")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("EBIN (EORTC-1612-MG)")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .evidenceUrls(Lists.newArrayList("https://www.google.com/#q=FDA"))
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.ICLUSION),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null))
                .build());

        trialsOnLabel.add(onLabelBuilder.event("High tumor mutation load")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("BASKET OF BASKETS (VHIO17002)")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .evidenceUrls(Lists.newArrayList("https://www.google.com/#q=FDA"))
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.ICLUSION), ProtectEvidenceType.SIGNATURE, null))
                .build());

        trialsOnLabel.add(onLabelBuilder.event("High tumor mutation load")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("DRUP")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .evidenceUrls(Lists.newArrayList("https://www.google.com/#q=FDA"))
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.ICLUSION), ProtectEvidenceType.SIGNATURE, null))
                .build());

        trialsOnLabel.add(onLabelBuilder.event("High tumor mutation load")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("KEYNOTE-158")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .evidenceUrls(Lists.newArrayList("https://www.google.com/#q=FDA"))
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.ICLUSION), ProtectEvidenceType.SIGNATURE, null))
                .build());

        trialsOnLabel.add(onLabelBuilder.gene("PTEN")
                .transcript("123")
                .isCanonical(true)
                .eventIsHighDriver(true)
                .event("partial loss")
                .germline(false)
                .reported(true)
                .treatment("DRUP")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .evidenceUrls(Lists.newArrayList("https://www.google.com/#q=FDA"))
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.ICLUSION), ProtectEvidenceType.DELETION, null))
                .build());

        return trialsOnLabel;
    }

    @NotNull
    private static List<ProtectEvidence> createCOLO829OffLabelEvidence() {
        List<ProtectEvidence> evidenceItemsOffLabel = Lists.newArrayList();

        ImmutableProtectEvidence.Builder offLabelBuilder = ImmutableProtectEvidence.builder().onLabel(false);

        evidenceItemsOffLabel.add(offLabelBuilder.gene("BRAF")
                .transcript("123")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("Bevacizumab")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .evidenceUrls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/19603024",
                        "http://www.ncbi.nlm.nih.gov/pubmed/19571295"))
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.VICC_CGI),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null))
                .build());

        evidenceItemsOffLabel.add(offLabelBuilder.gene("BRAF")
                .transcript("123")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("CI-1040")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.VICC_CGI),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null))
                .evidenceUrls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/18682506",
                        "http://www.ncbi.nlm.nih.gov/pubmed/21882184"))
                .build());

        evidenceItemsOffLabel.add(offLabelBuilder.gene("BRAF")
                .transcript("123")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("Cetuximab")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.VICC_CGI),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null))
                .evidenceUrls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/19001320",
                        "http://www.ncbi.nlm.nih.gov/pubmed/20619739",
                        "http://www.ncbi.nlm.nih.gov/pubmed/19884556",
                        "http://www.ncbi.nlm.nih.gov/pubmed/21163703",
                        "http://www.ncbi.nlm.nih.gov/pubmed/23325582",
                        "http://www.ncbi.nlm.nih.gov/pubmed/25666295"))
                .build());

        evidenceItemsOffLabel.add(offLabelBuilder.gene("BRAF")
                .transcript("123")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("Cetuximab + Irinotecan + Vemurafenib")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.VICC_CIVIC),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null))
                .evidenceUrls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/27729313"))
                .build());

        evidenceItemsOffLabel.add(offLabelBuilder.gene("BRAF")
                .transcript("123")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("Fluorouracil")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.VICC_CIVIC),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null))
                .evidenceUrls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/19603024"))
                .build());

        evidenceItemsOffLabel.add(offLabelBuilder.gene("BRAF")
                .transcript("123")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("Irinotecan")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.VICC_CIVIC),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null))
                .evidenceUrls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/19603024"))
                .build());

        evidenceItemsOffLabel.add(offLabelBuilder.gene("BRAF")
                .transcript("123")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("Oxaliplatin")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.VICC_CIVIC),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null))
                .evidenceUrls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/19603024"))
                .build());

        evidenceItemsOffLabel.add(offLabelBuilder.gene("BRAF")
                .transcript("123")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("Panitumumab")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.VICC_CGI, Knowledgebase.VICC_CIVIC),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null))
                .evidenceUrls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/21163703",
                        "http://www.ncbi.nlm.nih.gov/pubmed/23325582",
                        "http://www.ncbi.nlm.nih.gov/pubmed/19001320",
                        "http://www.ncbi.nlm.nih.gov/pubmed/20619739"))
                .build());

        evidenceItemsOffLabel.add(offLabelBuilder.gene("BRAF")
                .transcript("123")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("Selumetinib")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.VICC_CIVIC),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null))
                .evidenceUrls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/22492957"))
                .build());

        evidenceItemsOffLabel.add(offLabelBuilder.gene("BRAF")
                .transcript("123")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("Sorafenib")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.VICC_CIVIC),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null))
                .evidenceUrls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/18682506",
                        "http://www.ncbi.nlm.nih.gov/pubmed/21882184"))
                .build());

        evidenceItemsOffLabel.add(offLabelBuilder.gene("BRAF")
                .transcript("123")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("Vemurafenib")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.VICC_CGI), ProtectEvidenceType.CODON_MUTATION, 600))
                .evidenceUrls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/26287849"))
                .build());

        evidenceItemsOffLabel.add(offLabelBuilder.gene("PTEN")
                .transcript("123")
                .isCanonical(true)
                .event("partial loss")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("Anti-EGFR monoclonal antibody")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.VICC_CGI), ProtectEvidenceType.DELETION, null))
                .evidenceUrls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/21163703",
                        "http://www.ncbi.nlm.nih.gov/pubmed/19398573"))
                .build());

        evidenceItemsOffLabel.add(offLabelBuilder.gene("PTEN")
                .transcript("123")
                .isCanonical(true)
                .event("partial loss")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("Cetuximab")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.VICC_CIVIC), ProtectEvidenceType.DELETION, null))
                .evidenceUrls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/21163703"))
                .build());

        evidenceItemsOffLabel.add(offLabelBuilder.gene("PTEN")
                .transcript("123")
                .isCanonical(true)
                .event("partial loss")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("Everolimus")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.VICC_CIVIC), ProtectEvidenceType.DELETION, null))
                .evidenceUrls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/23989949"))
                .build());

        evidenceItemsOffLabel.add(offLabelBuilder.gene("PTEN")
                .transcript("123")
                .isCanonical(true)
                .event("partial loss")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("Lapatinib + Trastuzumab")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.VICC_CIVIC), ProtectEvidenceType.DELETION, null))
                .evidenceUrls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/25300346"))
                .build());

        evidenceItemsOffLabel.add(offLabelBuilder.gene("PTEN")
                .transcript("123")
                .isCanonical(true)
                .event("partial loss")
                .eventIsHighDriver(true)
                .germline(false)
                .reported(true)
                .treatment("Trastuzumab")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .protectSources(createTestProtectSource(Sets.newHashSet(Knowledgebase.VICC_CIVIC), ProtectEvidenceType.DELETION, null))
                .evidenceUrls(Lists.newArrayList("http://www.ncbi.nlm.nih.gov/pubmed/24387334",
                        "http://www.ncbi.nlm.nih.gov/pubmed/20813970"))
                .build());

        return evidenceItemsOffLabel;
    }

    @NotNull
    private static List<ReportableVariant> createCOLO829SomaticVariants(boolean forceCDKN2AVariantToBeGermline) {
        // TODO Fill in canonicalEffect
        ReportableVariant variant1 = ImmutableReportableVariant.builder()
                .source(ReportableVariantSource.SOMATIC)
                .gene("BRAF")
                .transcript("123")
                .isCanonical(true)
                .genotypeStatus(GenotypeStatus.HOM_REF)
                .chromosome("7")
                .position(140453136)
                .ref("T")
                .alt("A")
                .type(VariantType.SNP)
                .otherReportedEffects(Strings.EMPTY)
                .canonicalTranscript("ENST00000288602")
                .canonicalEffect(Strings.EMPTY)
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .canonicalHgvsCodingImpact("c.1799T>A")
                .canonicalHgvsProteinImpact("p.Val600Glu")
                .alleleReadCount(150)
                .totalReadCount(221)
                .alleleCopyNumber(4.09281)
                .totalCopyNumber(6.01)
                .minorAlleleCopyNumber(2.01)
                .hotspot(Hotspot.HOTSPOT)
                .driverLikelihood(1D)
                .clonalLikelihood(1D)
                .biallelic(false)
                .build();

        ReportableVariant variant2 = ImmutableReportableVariant.builder()
                .source(forceCDKN2AVariantToBeGermline ? ReportableVariantSource.GERMLINE : ReportableVariantSource.SOMATIC)
                .gene("CDKN2A")
                .transcript("123")
                .isCanonical(true)
                .genotypeStatus(GenotypeStatus.HOM_REF)
                .chromosome("9")
                .position(21971153)
                .ref("CCG")
                .alt("C")
                .type(VariantType.INDEL)
                .otherReportedEffects(Strings.EMPTY)
                .canonicalTranscript("ENST00000498124")
                .canonicalEffect(Strings.EMPTY)
                .canonicalCodingEffect(CodingEffect.NONSENSE_OR_FRAMESHIFT)
                .canonicalHgvsCodingImpact("c.203_204delCG")
                .canonicalHgvsProteinImpact("p.Ala68fs")
                .alleleReadCount(99)
                .totalReadCount(99)
                .alleleCopyNumber(2.0)
                .minorAlleleCopyNumber(0.0)
                .totalCopyNumber(2.0)
                .hotspot(Hotspot.NEAR_HOTSPOT)
                .clonalLikelihood(1D)
                .driverLikelihood(1D)
                .biallelic(true)
                .build();

        ReportableVariant variant3 = ImmutableReportableVariant.builder()
                .source(ReportableVariantSource.SOMATIC)
                .gene("TERT")
                .transcript("123")
                .isCanonical(true)
                .genotypeStatus(GenotypeStatus.HOM_REF)
                .chromosome("5")
                .position(1295228)
                .ref("GG")
                .alt("AA")
                .type(VariantType.MNP)
                .otherReportedEffects(Strings.EMPTY)
                .canonicalTranscript("ENST00000310581")
                .canonicalEffect("upstream_gene_variant")
                .canonicalCodingEffect(CodingEffect.NONE)
                .canonicalHgvsCodingImpact("c.-125_-124delCCinsTT")
                .canonicalHgvsProteinImpact(Strings.EMPTY)
                .alleleReadCount(56)
                .totalReadCount(65)
                .alleleCopyNumber(1.7404)
                .minorAlleleCopyNumber(0.0)
                .totalCopyNumber(2.0)
                .hotspot(Hotspot.HOTSPOT)
                .clonalLikelihood(1D)
                .driverLikelihood(1D)
                .biallelic(true)
                .build();

        ReportableVariant variant4 = ImmutableReportableVariant.builder()
                .source(ReportableVariantSource.SOMATIC)
                .gene("SF3B1")
                .transcript("123")
                .isCanonical(true)
                .genotypeStatus(GenotypeStatus.HOM_REF)
                .chromosome("2")
                .position(198266779)
                .ref("G")
                .alt("A")
                .type(VariantType.SNP)
                .otherReportedEffects(Strings.EMPTY)
                .canonicalTranscript("ENST00000335508")
                .canonicalEffect(Strings.EMPTY)
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .canonicalHgvsCodingImpact("c.2153C>T")
                .canonicalHgvsProteinImpact("p.Pro718Leu")
                .alleleReadCount(74)
                .totalReadCount(111)
                .alleleCopyNumber(2.026722)
                .minorAlleleCopyNumber(1.0)
                .totalCopyNumber(3.02)
                .hotspot(Hotspot.NON_HOTSPOT)
                .clonalLikelihood(1D)
                .driverLikelihood(0.1459)
                .biallelic(false)
                .build();

        ReportableVariant variant5 = ImmutableReportableVariant.builder()
                .source(ReportableVariantSource.SOMATIC)
                .gene("TP63")
                .transcript("123")
                .isCanonical(true)
                .genotypeStatus(GenotypeStatus.HOM_REF)
                .chromosome("3")
                .position(189604330)
                .ref("G")
                .alt("T")
                .type(VariantType.SNP)
                .otherReportedEffects(Strings.EMPTY)
                .canonicalTranscript("ENST00000264731")
                .canonicalEffect(Strings.EMPTY)
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .canonicalHgvsCodingImpact("c.1497G>T")
                .canonicalHgvsProteinImpact("p.Met499Ile")
                .alleleReadCount(47)
                .totalReadCount(112)
                .alleleCopyNumber(1.678764)
                .minorAlleleCopyNumber(1.97)
                .totalCopyNumber(3.98)
                .hotspot(Hotspot.NON_HOTSPOT)
                .clonalLikelihood(1D)
                .driverLikelihood(0)
                .biallelic(false)
                .build();

        return Lists.newArrayList(variant1, variant2, variant3, variant4, variant5);
    }

    @NotNull
    private static List<ReportableGainLoss> createCOLO829GainsLosses() {
        ReportableGainLoss gainLoss1 = ImmutableReportableGainLoss.builder()
                .chromosome("10")
                .chromosomeBand("q23.31")
                .gene("PTEN")
                .transcript("123")
                .isCanonical(true)
                .minCopies(0)
                .maxCopies(2)
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
                .reportedType(KnownFusionType.KNOWN_PAIR.toString())
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
                .reportedType(KnownFusionType.PROMISCUOUS_5.toString())
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
        ReportableGeneDisruption disruption1 = ImmutableReportableGeneDisruption.builder()
                .location("10q23.31")
                .gene("PTEN")
                .range("Intron 5 -> Intron 6")
                .type("DEL")
                .junctionCopyNumber(2.012)
                .undisruptedCopyNumber(0.0)
                .firstAffectedExon(5)
                .svId(1)
                .clusterId(2)
                .build();

        return Lists.newArrayList(disruption1);
    }

    @NotNull
    private static List<ReportableHomozygousDisruption> createTestHomozygousDisruptions() {
        return Lists.newArrayList(ImmutableReportableHomozygousDisruption.builder()
                .chromosome("8")
                .chromosomeBand("p22")
                .gene("SGCZ")
                .transcript("123")
                .isCanonical(true)
                .build());
    }

    @NotNull
    private static List<PeachGenotype> createTestPeachGenotypes() {
        return Lists.newArrayList(ImmutablePeachGenotype.builder()
                .gene("DPYD")
                .haplotype("*1_HOM")
                .function("Normal Function")
                .linkedDrugs("5-Fluorouracil;Capecitabine;Tegafur")
                .urlPrescriptionInfo("https://www.pharmgkb.org/chemical/PA128406956/guidelineAnnotation/PA166104939;"
                        + "https://www.pharmgkb.org/chemical/PA448771/guidelineAnnotation/PA166104963;"
                        + "https://www.pharmgkb.org/chemical/PA452620/guidelineAnnotation/PA166104944")
                .panelVersion("PGx_min_DPYD_v1.2")
                .repoVersion("1.4")
                .build());
    }

    @NotNull
    private static List<AnnotatedVirus> createTestAnnotatedViruses() {
        return Lists.newArrayList(VirusTestFactory.testAnnotatedVirusBuilder()
                .name("Human papillomavirus type 16")
                .integrations(2)
                .interpretation("HPV")
                .reported(true)
                .build());
    }
}