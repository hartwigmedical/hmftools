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
import com.hartwig.hmftools.patientreporter.cfreport.MathUtil;
import com.hartwig.hmftools.patientreporter.cfreport.data.TumorPurity;

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
        String pipelineVersion = "5.28";
        double averageTumorPloidy = 3.1;
        int tumorMutationalLoad = 186;
        double tumorMutationalBurden = 13.7205;
        double microsatelliteIndelsPerMb = 0.1172;
        double chordHrdValue = 0D;
        ChordStatus chordStatus = ChordStatus.HR_PROFICIENT;
        String reportDate = DataUtil.formatDate(LocalDate.now());
        double impliedPurityPercentage = MathUtil.mapPercentage(config.impliedTumorPurity(), TumorPurity.RANGE_MIN, TumorPurity.RANGE_MAX);

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

        String summaryWithoutGermlineLowPurity = "Melanoma sample showing:\n"
                + " - activating BRAF mutation that is associated with response to BRAF-inhibitors (in combination with a MEK-inhibitor)\n"
                + " - complete inactivation of CDKN2A, indicating potential benefit of CDK4/6 inhibitors\n"
                + " - complete inactivation/loss of PTEN likely resulting in an activation of the PI3K-AKT-mTOR pathway "
                + "and indicating potential benefit of mTOR/PI3K inhibitors\n"
                + " - high mutational burden (mutational load (ML) of 180, tumor mutation burden (TMB) of 13.6) that is "
                + "potentially associated with an increased response rate to checkpoint inhibitor immunotherapy\n"
                + "Due to the lower tumor purity (" + DataUtil.formatPercentage(impliedPurityPercentage) + ") potential (subclonal) "
                + "DNA aberrations might not have been detected using this test. This result should therefore be considered with caution.";

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
            if (config.reportGermline()) {
                clinicalSummary = summaryWithGermline;
            } else {
                if (purpleQCStatus == PurpleQCStatus.FAIL_NO_TUMOR || purpleQCStatus == PurpleQCStatus.WARN_LOW_PURITY) {
                    clinicalSummary = summaryWithoutGermlineLowPurity;
                } else {
                    clinicalSummary = summaryWithoutGermline;
                }
            }
        } else {
            clinicalSummary = Strings.EMPTY;
        }

        String specialremark = "This is a special remark";

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
                .specialRemark(specialremark)
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
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("1"), ChromosomeArm.P_ARM, 2.5772792739090105));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("1"), ChromosomeArm.Q_ARM, 3.923536090197291));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("2"), ChromosomeArm.P_ARM, 3.016271458549662));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("2"), ChromosomeArm.Q_ARM, 3.02120002022585));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("3"), ChromosomeArm.P_ARM, 3.5910709809245502));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("3"), ChromosomeArm.Q_ARM, 4.000879324279213));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("4"), ChromosomeArm.P_ARM, 2.0210999604946176));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("4"), ChromosomeArm.Q_ARM, 3.8453841055332885));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("5"), ChromosomeArm.P_ARM, 2.0000481443928493));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("5"), ChromosomeArm.Q_ARM, 2.0096989504909413));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("6"), ChromosomeArm.P_ARM, 3.8479440436530545));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("6"), ChromosomeArm.Q_ARM, 2.913192184870031));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("7"), ChromosomeArm.P_ARM, 4.024620078272729));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("7"), ChromosomeArm.Q_ARM, 4.171334679722508));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("8"), ChromosomeArm.P_ARM, 3.3325999264957695));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("8"), ChromosomeArm.Q_ARM, 3.3429529791412795));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("9"), ChromosomeArm.P_ARM, 2.7291755064569365));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("9"), ChromosomeArm.Q_ARM, 3.6992000400581495));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("10"), ChromosomeArm.P_ARM, 2.5009802424099066));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("10"), ChromosomeArm.Q_ARM, 2.007094743290902));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("11"), ChromosomeArm.P_ARM, 3.1661435478205004));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("11"), ChromosomeArm.Q_ARM, 2.9098638260285616));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("12"), ChromosomeArm.P_ARM, 3.0115999171651855));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("12"), ChromosomeArm.Q_ARM, 3.003116732830778));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("13"), ChromosomeArm.P_ARM, 3.1564998196285714));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("13"), ChromosomeArm.Q_ARM, 3.146714774779385));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("14"), ChromosomeArm.P_ARM, 3.0139998277714284));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("14"), ChromosomeArm.Q_ARM, 3.0138455205847463));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("15"), ChromosomeArm.P_ARM, 3.7023997998702702));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("15"), ChromosomeArm.Q_ARM, 2.546611471333237));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("16"), ChromosomeArm.P_ARM, 3.1919712830460782));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("16"), ChromosomeArm.Q_ARM, 1.989521251230779));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("17"), ChromosomeArm.P_ARM, 2.9938998740100473));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("17"), ChromosomeArm.Q_ARM, 3.0477000530660465));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("18"), ChromosomeArm.P_ARM, 2.370931614063123));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("18"), ChromosomeArm.Q_ARM, 2.8490432529511334));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("19"), ChromosomeArm.P_ARM, 2.889021931658433));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("19"), ChromosomeArm.Q_ARM, 2.934100089054606));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("20"), ChromosomeArm.P_ARM, 4.013880905298536));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("20"), ChromosomeArm.Q_ARM, 4.008612196597953));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("21"), ChromosomeArm.P_ARM, 2.991999766033014));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("21"), ChromosomeArm.Q_ARM, 2.9982181161009325));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("22"), ChromosomeArm.P_ARM, 3.9915997247172412));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("22"), ChromosomeArm.Q_ARM, 3.983474946385728));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("X"), ChromosomeArm.P_ARM, 1.9496150872949336));
        cnPerChromosomeArm.add(buildCnPerChromosomeArmData(HumanChromosome.fromString("X"), ChromosomeArm.Q_ARM, 1.9547000205458256));
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
    private static ProtectSource createTestProtectSource(@NotNull Knowledgebase source, @NotNull String soureEvent,
            @NotNull Set<String> sourceUrls, @NotNull ProtectEvidenceType protectEvidenceType, @Nullable Integer rank,
            @NotNull Set<String> evidenceurls) {

        return ImmutableProtectSource.builder()
                .source(source)
                .sourceEvent(soureEvent)
                .sourceUrls(sourceUrls)
                .evidenceType(protectEvidenceType)
                .rangeRank(rank)
                .evidenceUrls(evidenceurls)
                .build();
    }

    @NotNull
    private static List<ProtectEvidence> createCOLO829TumorSpecificEvidence() {
        List<ProtectEvidence> evidenceItemsOnLabel = Lists.newArrayList();

        ImmutableProtectEvidence.Builder onLabelBuilder = ImmutableProtectEvidence.builder().onLabel(true).reported(true);
        evidenceItemsOnLabel.add(onLabelBuilder.gene("BRAF")
                .transcript("ENST00000288602")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .treatment("Cobimetinib + Vemurafenib")
                .level(EvidenceLevel.A)
                .direction(EvidenceDirection.RESPONSIVE)
                .protectSources(Lists.newArrayList(createTestProtectSource(Knowledgebase.VICC_CGI,
                        "BRAF BRAF:V600E",
                        Sets.newHashSet(),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null,
                        Sets.newHashSet("ttps://www.google.com/#q=FDA"))))
                .build());

        evidenceItemsOnLabel.add(onLabelBuilder.gene("BRAF")
                .transcript("ENST00000288602")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .treatment("Dabrafenib")
                .level(EvidenceLevel.A)
                .direction(EvidenceDirection.RESPONSIVE)
                .protectSources(Lists.newArrayList(createTestProtectSource(Knowledgebase.VICC_CGI,
                        "BRAF BRAF:V600E",
                        Sets.newHashSet(),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null,
                        Sets.newHashSet("https://www.google.com/#q=FDA"))))
                .build());

        evidenceItemsOnLabel.add(onLabelBuilder.gene("BRAF")
                .transcript("ENST00000288602")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .treatment("Dabrafenib + Trametinib")
                .level(EvidenceLevel.A)
                .direction(EvidenceDirection.RESPONSIVE)
                .protectSources(Lists.newArrayList(createTestProtectSource(Knowledgebase.VICC_CIVIC,
                                "BRAF V600E",
                                Sets.newHashSet(),
                                ProtectEvidenceType.HOTSPOT_MUTATION,
                                null,
                                Sets.newHashSet("http://www.ncbi.nlm.nih.gov/pubmed/25399551")),
                        createTestProtectSource(Knowledgebase.VICC_CGI,
                                "BRAF BRAF:V600E",
                                Sets.newHashSet(),
                                ProtectEvidenceType.HOTSPOT_MUTATION,
                                null,
                                Sets.newHashSet("https://www.google.com/#q=FDA"))))
                .build());

        evidenceItemsOnLabel.add(onLabelBuilder.gene("BRAF")
                .transcript("ENST00000288602")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .treatment("Trametinib")
                .level(EvidenceLevel.A)
                .direction(EvidenceDirection.RESPONSIVE)
                .protectSources(Lists.newArrayList(createTestProtectSource(Knowledgebase.VICC_CGI,
                        "BRAF BRAF:V600E",
                        Sets.newHashSet(),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null,
                        Sets.newHashSet("https://www.google.com/#q=FDA"))))
                .build());

        evidenceItemsOnLabel.add(onLabelBuilder.gene("BRAF")
                .transcript("ENST00000288602")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .treatment("Vemurafenib")
                .level(EvidenceLevel.A)
                .direction(EvidenceDirection.RESPONSIVE)
                .protectSources(Lists.newArrayList(createTestProtectSource(Knowledgebase.VICC_CGI,
                                "BRAF BRAF:V600G",
                                Sets.newHashSet(),
                                ProtectEvidenceType.HOTSPOT_MUTATION,
                                null,
                                Sets.newHashSet("https://www.google.com/#q=NCCN")),
                        createTestProtectSource(Knowledgebase.VICC_CGI,
                                "BRAF BRAF:V600E",
                                Sets.newHashSet(),
                                ProtectEvidenceType.HOTSPOT_MUTATION,
                                null,
                                Sets.newHashSet("https://www.google.com/#q=FDA")),
                        createTestProtectSource(Knowledgebase.VICC_CIVIC,
                                "BRAF V600E",
                                Sets.newHashSet(),
                                ProtectEvidenceType.HOTSPOT_MUTATION,
                                null,
                                Sets.newHashSet("http://www.ncbi.nlm.nih.gov/pubmed/21639808"))))
                .build());

        evidenceItemsOnLabel.add(onLabelBuilder.gene("BRAF")
                .transcript("ENST00000288602")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .treatment("RO4987655")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .protectSources(Lists.newArrayList(createTestProtectSource(Knowledgebase.VICC_CIVIC,
                        "BRAF V600",
                        Sets.newHashSet(),
                        ProtectEvidenceType.CODON_MUTATION,
                        600,
                        Sets.newHashSet("http://www.ncbi.nlm.nih.gov/pubmed/24947927"))))
                .build());

        evidenceItemsOnLabel.add(onLabelBuilder.gene("PTEN")
                .transcript("ENST00000371953")
                .isCanonical(true)
                .event("partial loss")
                .eventIsHighDriver(true)
                .germline(false)
                .treatment("Buparlisib + Carboplatin + Paclitaxel")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .protectSources(Lists.newArrayList(createTestProtectSource(Knowledgebase.VICC_CIVIC,
                        "PTEN LOSS",
                        Sets.newHashSet(),
                        ProtectEvidenceType.DELETION,
                        null,
                        Sets.newHashSet("http://www.ncbi.nlm.nih.gov/pubmed/25672916"))))
                .build());
        return evidenceItemsOnLabel;
    }

    @NotNull
    private static List<ProtectEvidence> createCOLO829ClinicalTrials() {
        List<ProtectEvidence> trialsOnLabel = Lists.newArrayList();
        ImmutableProtectEvidence.Builder trialBuilder = ImmutableProtectEvidence.builder().onLabel(true).reported(true);

        trialsOnLabel.add(trialBuilder.gene(null)
                .transcript(null)
                .isCanonical(false)
                .event("High tumor mutation load")
                .eventIsHighDriver(null)
                .germline(false)
                .treatment("BASKET OF BASKETS (VHIO17002)")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .protectSources(Lists.newArrayList(createTestProtectSource(Knowledgebase.ICLUSION,
                        "TumMutLoad_HIGH",
                        Sets.newHashSet("https://trial-eye.com/hmf/11087"),
                        ProtectEvidenceType.SIGNATURE,
                        null,
                        Sets.newHashSet())))
                .build());

        trialsOnLabel.add(trialBuilder.gene(null)
                .transcript(null)
                .isCanonical(false)
                .event("High tumor mutation load")
                .eventIsHighDriver(null)
                .germline(false)
                .treatment("DRUP")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .protectSources(Lists.newArrayList(createTestProtectSource(Knowledgebase.ICLUSION,
                        "TumMutLoad_HIGH",
                        Sets.newHashSet("https://trial-eye.com/hmf/10299"),
                        ProtectEvidenceType.SIGNATURE,
                        null,
                        Sets.newHashSet())))
                .build());

        trialsOnLabel.add(trialBuilder.gene(null)
                .transcript(null)
                .isCanonical(false)
                .event("High tumor mutation load")
                .eventIsHighDriver(null)
                .germline(false)
                .treatment("KEYNOTE-158")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .protectSources(Lists.newArrayList(createTestProtectSource(Knowledgebase.ICLUSION,
                        "TumMutLoad_HIGH",
                        Sets.newHashSet("https://trial-eye.com/hmf/4866"),
                        ProtectEvidenceType.SIGNATURE,
                        null,
                        Sets.newHashSet())))
                .build());

        trialsOnLabel.add(trialBuilder.gene("BRAF")
                .transcript("ENST00000288602")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .treatment("Array 818-103")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .protectSources(Lists.newArrayList(createTestProtectSource(Knowledgebase.ICLUSION,
                        "BRAF V600",
                        Sets.newHashSet("https://trial-eye.com/hmf/13054"),
                        ProtectEvidenceType.CODON_MUTATION,
                        600,
                        Sets.newHashSet())))
                .build());

        trialsOnLabel.add(trialBuilder.gene("BRAF")
                .transcript("ENST00000288602")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .treatment("DRUP")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .protectSources(Lists.newArrayList(createTestProtectSource(Knowledgebase.ICLUSION,
                                "BRAF ACTIVATING MUTATION",
                                Sets.newHashSet("https://trial-eye.com/hmf/10299"),
                                ProtectEvidenceType.ACTIVATION,
                                null,
                                Sets.newHashSet()),
                        createTestProtectSource(Knowledgebase.ICLUSION,
                                "BRAF V600",
                                Sets.newHashSet("https://trial-eye.com/hmf/10299"),
                                ProtectEvidenceType.CODON_MUTATION,
                                600,
                                Sets.newHashSet())))
                .build());

        trialsOnLabel.add(trialBuilder.gene("BRAF")
                .transcript("ENST00000288602")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .treatment("EBIN (EORTC-1612-MG)")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .protectSources(Lists.newArrayList(createTestProtectSource(Knowledgebase.ICLUSION,
                        "BRAF V600",
                        Sets.newHashSet("https://trial-eye.com/hmf/11284"),
                        ProtectEvidenceType.CODON_MUTATION,
                        600,
                        Sets.newHashSet())))
                .build());

        trialsOnLabel.add(trialBuilder.gene("BRAF")
                .transcript("ENST00000288602")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .treatment("NASAM")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .protectSources(Lists.newArrayList(createTestProtectSource(Knowledgebase.ICLUSION,
                        "BRAF V600E",
                        Sets.newHashSet("https://trial-eye.com/hmf/14995"),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null,
                        Sets.newHashSet())))
                .build());

        trialsOnLabel.add(trialBuilder.gene("PTEN")
                .transcript("ENST00000371953")
                .isCanonical(true)
                .event("partial loss")
                .eventIsHighDriver(true)
                .germline(false)
                .treatment("DRUP")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .protectSources(Lists.newArrayList(createTestProtectSource(Knowledgebase.ICLUSION,
                                "PTEN LOSS",
                                Sets.newHashSet("https://trial-eye.com/hmf/10299\""),
                                ProtectEvidenceType.DELETION,
                                null,
                                Sets.newHashSet()),
                        createTestProtectSource(Knowledgebase.ICLUSION,
                                "PTEN INACTIVATION MUTATION",
                                Sets.newHashSet("https://trial-eye.com/hmf/10299"),
                                ProtectEvidenceType.INACTIVATION,
                                null,
                                Sets.newHashSet())))
                .build());
        return trialsOnLabel;
    }

    @NotNull
    private static List<ProtectEvidence> createCOLO829OffLabelEvidence() {
        List<ProtectEvidence> evidenceItemsOffLabel = Lists.newArrayList();

        ImmutableProtectEvidence.Builder offLabelBuilder = ImmutableProtectEvidence.builder().onLabel(false).reported(true);

        evidenceItemsOffLabel.add(offLabelBuilder.gene("BRAF")
                .transcript("ENST00000288602")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .treatment("Bevacizumab")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .protectSources(Lists.newArrayList(createTestProtectSource(Knowledgebase.VICC_CIVIC,
                        "BRAF V600E",
                        Sets.newHashSet(),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null,
                        Sets.newHashSet("http://www.ncbi.nlm.nih.gov/pubmed/19571295", "http://www.ncbi.nlm.nih.gov/pubmed/19603024"))))
                .build());

        evidenceItemsOffLabel.add(offLabelBuilder.gene("BRAF")
                .transcript("ENST00000288602")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .treatment("CI-1040")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .protectSources(Lists.newArrayList(createTestProtectSource(Knowledgebase.VICC_CIVIC,
                        "BRAF V600E",
                        Sets.newHashSet(),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null,
                        Sets.newHashSet("http://www.ncbi.nlm.nih.gov/pubmed/21882184", "http://www.ncbi.nlm.nih.gov/pubmed/18682506"))))
                .build());

        evidenceItemsOffLabel.add(offLabelBuilder.gene("BRAF")
                .transcript("ENST00000288602")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .treatment("Cetuximab")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .protectSources(Lists.newArrayList(createTestProtectSource(Knowledgebase.VICC_CGI,
                                "BRAF BRAF:V600E",
                                Sets.newHashSet(),
                                ProtectEvidenceType.HOTSPOT_MUTATION,
                                null,
                                Sets.newHashSet("http://www.ncbi.nlm.nih.gov/pubmed/20619739",
                                        "http://www.ncbi.nlm.nih.gov/pubmed/21163703",
                                        "http://www.ncbi.nlm.nih.gov/pubmed/23325582")),
                        createTestProtectSource(Knowledgebase.VICC_CIVIC,
                                "BRAF V600E",
                                Sets.newHashSet(),
                                ProtectEvidenceType.HOTSPOT_MUTATION,
                                null,
                                Sets.newHashSet("http://www.ncbi.nlm.nih.gov/pubmed/19001320",
                                        "http://www.ncbi.nlm.nih.gov/pubmed/20619739",
                                        "http://www.ncbi.nlm.nih.gov/pubmed/19884556",
                                        "http://www.ncbi.nlm.nih.gov/pubmed/25666295"))))
                .build());

        evidenceItemsOffLabel.add(offLabelBuilder.gene("BRAF")
                .transcript("ENST00000288602")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .treatment("Cetuximab + Irinotecan + Vemurafenib")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .protectSources(Lists.newArrayList(createTestProtectSource(Knowledgebase.VICC_CIVIC,
                        "BRAF V600E",
                        Sets.newHashSet(),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null,
                        Sets.newHashSet("http://www.ncbi.nlm.nih.gov/pubmed/27729313"))))
                .build());

        evidenceItemsOffLabel.add(offLabelBuilder.gene("BRAF")
                .transcript("ENST00000288602")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .treatment("Fluorouracil")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .protectSources(Lists.newArrayList(createTestProtectSource(Knowledgebase.VICC_CIVIC,
                        "BRAF V600E",
                        Sets.newHashSet(),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null,
                        Sets.newHashSet("http://www.ncbi.nlm.nih.gov/pubmed/19603024"))))
                .build());

        evidenceItemsOffLabel.add(offLabelBuilder.gene("BRAF")
                .transcript("ENST00000288602")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .treatment("Irinotecan")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .protectSources(Lists.newArrayList(createTestProtectSource(Knowledgebase.VICC_CIVIC,
                        "BRAF V600E",
                        Sets.newHashSet(),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null,
                        Sets.newHashSet("http://www.ncbi.nlm.nih.gov/pubmed/19603024"))))
                .build());

        evidenceItemsOffLabel.add(offLabelBuilder.gene("BRAF")
                .transcript("ENST00000288602")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .treatment("Oxaliplatin")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .protectSources(Lists.newArrayList(createTestProtectSource(Knowledgebase.VICC_CIVIC,
                        "BRAF V600E",
                        Sets.newHashSet(),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null,
                        Sets.newHashSet("http://www.ncbi.nlm.nih.gov/pubmed/19603024"))))
                .build());

        evidenceItemsOffLabel.add(offLabelBuilder.gene("BRAF")
                .transcript("ENST00000288602")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .treatment("Panitumumab")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .protectSources(Lists.newArrayList(createTestProtectSource(Knowledgebase.VICC_CGI,
                                "BRAF BRAF:V600E",
                                Sets.newHashSet(),
                                ProtectEvidenceType.HOTSPOT_MUTATION,
                                null,
                                Sets.newHashSet("http://www.ncbi.nlm.nih.gov/pubmed/20619739",
                                        "http://www.ncbi.nlm.nih.gov/pubmed/21163703",
                                        "http://www.ncbi.nlm.nih.gov/pubmed/23325582")),
                        createTestProtectSource(Knowledgebase.VICC_CIVIC,
                                "BRAF V600",
                                Sets.newHashSet(),
                                ProtectEvidenceType.CODON_MUTATION,
                                600,
                                Sets.newHashSet("http://www.ncbi.nlm.nih.gov/pubmed/23325582")),
                        createTestProtectSource(Knowledgebase.VICC_CIVIC,
                                "BRAF V600E",
                                Sets.newHashSet(),
                                ProtectEvidenceType.HOTSPOT_MUTATION,
                                null,
                                Sets.newHashSet("http://www.ncbi.nlm.nih.gov/pubmed/19001320"))))
                .build());

        evidenceItemsOffLabel.add(offLabelBuilder.gene("BRAF")
                .transcript("ENST00000288602")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .treatment("Selumetinib")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .protectSources(Lists.newArrayList(createTestProtectSource(Knowledgebase.VICC_CIVIC,
                        "BRAF V600E",
                        Sets.newHashSet(),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null,
                        Sets.newHashSet("http://www.ncbi.nlm.nih.gov/pubmed/22492957"))))
                .build());

        evidenceItemsOffLabel.add(offLabelBuilder.gene("BRAF")
                .transcript("ENST00000288602")
                .isCanonical(true)
                .event("p.Val600Glu")
                .eventIsHighDriver(true)
                .germline(false)
                .treatment("Sorafenib")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .protectSources(Lists.newArrayList(createTestProtectSource(Knowledgebase.VICC_CIVIC,
                        "BRAF V600E",
                        Sets.newHashSet(),
                        ProtectEvidenceType.HOTSPOT_MUTATION,
                        null,
                        Sets.newHashSet("http://www.ncbi.nlm.nih.gov/pubmed/21882184", "http://www.ncbi.nlm.nih.gov/pubmed/18682506"))))
                .build());

        evidenceItemsOffLabel.add(offLabelBuilder.gene("PTEN")
                .transcript("ENST00000371953")
                .isCanonical(true)
                .event("partial loss")
                .eventIsHighDriver(true)
                .germline(false)
                .treatment("Anti-EGFR monoclonal antibody")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .protectSources(Lists.newArrayList(createTestProtectSource(Knowledgebase.VICC_CGI,
                                "PTEN PTEN oncogenic mutation",
                                Sets.newHashSet(),
                                ProtectEvidenceType.INACTIVATION,
                                null,
                                Sets.newHashSet("http://www.ncbi.nlm.nih.gov/pubmed/21163703", "http://www.ncbi.nlm.nih.gov/pubmed/19398573")),
                        createTestProtectSource(Knowledgebase.VICC_CGI,
                                "PTEN PTEN deletion",
                                Sets.newHashSet(),
                                ProtectEvidenceType.DELETION,
                                null,
                                Sets.newHashSet("http://www.ncbi.nlm.nih.gov/pubmed/21163703",
                                        "http://www.ncbi.nlm.nih.gov/pubmed/19398573"))))
                .build());

        evidenceItemsOffLabel.add(offLabelBuilder.gene("PTEN")
                .transcript("ENST00000371953")
                .isCanonical(true)
                .event("partial loss")
                .eventIsHighDriver(true)
                .germline(false)
                .treatment("Cetuximab")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .protectSources(Lists.newArrayList(createTestProtectSource(Knowledgebase.VICC_CIVIC,
                        "PTEN LOSS",
                        Sets.newHashSet(),
                        ProtectEvidenceType.DELETION,
                        null,
                        Sets.newHashSet("http://www.ncbi.nlm.nih.gov/pubmed/21163703"))))
                .build());

        evidenceItemsOffLabel.add(offLabelBuilder.gene("PTEN")
                .transcript("ENST00000371953")
                .isCanonical(true)
                .event("partial loss")
                .eventIsHighDriver(true)
                .germline(false)
                .treatment("Everolimus")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .protectSources(Lists.newArrayList(createTestProtectSource(Knowledgebase.VICC_CIVIC,
                        "PTEN LOSS",
                        Sets.newHashSet(),
                        ProtectEvidenceType.DELETION,
                        null,
                        Sets.newHashSet("http://www.ncbi.nlm.nih.gov/pubmed/23989949"))))
                .build());

        evidenceItemsOffLabel.add(offLabelBuilder.gene("PTEN")
                .transcript("ENST00000371953")
                .isCanonical(true)
                .event("partial loss")
                .eventIsHighDriver(true)
                .germline(false)
                .treatment("Lapatinib + Trastuzumab")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .protectSources(Lists.newArrayList(createTestProtectSource(Knowledgebase.VICC_CIVIC,
                        "PTEN LOSS",
                        Sets.newHashSet(),
                        ProtectEvidenceType.DELETION,
                        null,
                        Sets.newHashSet("http://www.ncbi.nlm.nih.gov/pubmed/25300346"))))
                .build());

        evidenceItemsOffLabel.add(offLabelBuilder.gene("PTEN")
                .transcript("ENST00000371953")
                .isCanonical(true)
                .event("partial loss")
                .eventIsHighDriver(true)
                .germline(false)
                .treatment("Trastuzumab")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESISTANT)
                .protectSources(Lists.newArrayList(createTestProtectSource(Knowledgebase.VICC_CIVIC,
                                "PTEN LOSS",
                                Sets.newHashSet(),
                                ProtectEvidenceType.DELETION,
                                null,
                                Sets.newHashSet("http://www.ncbi.nlm.nih.gov/pubmed/20813970")),
                        createTestProtectSource(Knowledgebase.VICC_CIVIC,
                                "PTEN LOSS",
                                Sets.newHashSet(),
                                ProtectEvidenceType.DELETION,
                                null,
                                Sets.newHashSet("http://www.ncbi.nlm.nih.gov/pubmed/24387334"))))
                .build());
        return evidenceItemsOffLabel;
    }

    @NotNull
    private static List<ReportableVariant> createCOLO829SomaticVariants(boolean forceCDKN2AVariantToBeGermline) {
        ReportableVariant variant1 = ImmutableReportableVariant.builder()
                .source(ReportableVariantSource.SOMATIC)
                .gene("BRAF")
                .transcript("ENST00000288602")
                .isCanonical(true)
                .genotypeStatus(GenotypeStatus.HOM_REF)
                .chromosome("7")
                .position(140453136)
                .ref("A")
                .alt("T")
                .type(VariantType.SNP)
                .otherReportedEffects(Strings.EMPTY)
                .canonicalTranscript("ENST00000288602")
                .canonicalEffect("missense_variant")
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .canonicalHgvsCodingImpact("c.1799T>A")
                .canonicalHgvsProteinImpact("p.Val600Glu")
                .alleleReadCount(150)
                .totalReadCount(221)
                .alleleCopyNumber(4.09962)
                .totalCopyNumber(6.02)
                .minorAlleleCopyNumber(2.01)
                .hotspot(Hotspot.HOTSPOT)
                .driverLikelihood(1D)
                .clonalLikelihood(1D)
                .biallelic(false)
                .build();

        ReportableVariant variant2 = ImmutableReportableVariant.builder()
                .source(forceCDKN2AVariantToBeGermline ? ReportableVariantSource.GERMLINE : ReportableVariantSource.SOMATIC)
                .gene("CDKN2A (p16)")
                .transcript("ENST00000498124")
                .isCanonical(true)
                .genotypeStatus(GenotypeStatus.HOM_REF)
                .chromosome("9")
                .position(21971153)
                .ref("CCG")
                .alt("C")
                .type(VariantType.INDEL)
                .otherReportedEffects("ENST00000579755|c.246_247delCG|p.Gly83fs|frameshift_variant|NONSENSE_OR_FRAMESHIFT")
                .canonicalTranscript("ENST00000498124")
                .canonicalEffect("frameshift_variant")
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
                .transcript("ENST00000310581")
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
                .canonicalHgvsCodingImpact(Strings.EMPTY)
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
                .localPhaseSet(4410)
                .build();

        ReportableVariant variant4 = ImmutableReportableVariant.builder()
                .source(ReportableVariantSource.SOMATIC)
                .gene("SF3B1")
                .transcript("ENST00000335508")
                .isCanonical(true)
                .genotypeStatus(GenotypeStatus.HOM_REF)
                .chromosome("2")
                .position(198266779)
                .ref("G")
                .alt("A")
                .type(VariantType.SNP)
                .otherReportedEffects(Strings.EMPTY)
                .canonicalTranscript("ENST00000335508")
                .canonicalEffect("missense_variant")
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
                .driverLikelihood(0.1458)
                .biallelic(false)
                .build();

        ReportableVariant variant5 = ImmutableReportableVariant.builder()
                .source(ReportableVariantSource.SOMATIC)
                .gene("TP63")
                .transcript("ENST00000264731")
                .isCanonical(true)
                .genotypeStatus(GenotypeStatus.HOM_REF)
                .chromosome("3")
                .position(189604330)
                .ref("G")
                .alt("T")
                .type(VariantType.SNP)
                .otherReportedEffects(Strings.EMPTY)
                .canonicalTranscript("ENST00000264731")
                .canonicalEffect("missense_variant")
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
                .transcript("ENST00000371953")
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
                .junctionCopyNumber(2.0076)
                .undisruptedCopyNumber(0.0)
                .firstAffectedExon(5)
                .svId(121)
                .clusterId(67)
                .transcriptId("ENST00000371953")
                .isCanonical(true)
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
                .urlPrescriptionInfo("https://www.pharmgkb.org/chemical/PA128406956/guidelineAnnotation/PA166104939"
                        + "https://www.pharmgkb.org/chemical/PA448771/guidelineAnnotation/PA166104963"
                        + "https://www.pharmgkb.org/chemical/PA452620/guidelineAnnotation/PA166104944")
                .panelVersion("PGx_min_DPYD_v1.2")
                .repoVersion("1.6")
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