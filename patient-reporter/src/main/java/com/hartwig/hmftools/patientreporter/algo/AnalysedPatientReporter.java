package com.hartwig.hmftools.patientreporter.algo;

import java.io.File;
import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;
import java.util.Map;
import java.util.Optional;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.clinical.PatientPrimaryTumor;
import com.hartwig.hmftools.common.clinical.PatientPrimaryTumorFunctions;
import com.hartwig.hmftools.common.cuppa.CuppaDataFile;
import com.hartwig.hmftools.common.cuppa.interpretation.CuppaPrediction;
import com.hartwig.hmftools.common.cuppa.interpretation.CuppaPredictionFactory;
import com.hartwig.hmftools.common.cuppa.interpretation.CuppaReporting;
import com.hartwig.hmftools.common.cuppa.interpretation.CuppaReportingFactory;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.common.peach.PeachGenotypeFile;
import com.hartwig.hmftools.common.pipeline.PipelineVersionFile;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.rose.RoseConclusionFile;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariantSource;
import com.hartwig.hmftools.patientreporter.PatientReporterConfig;
import com.hartwig.hmftools.patientreporter.QsFormNumber;
import com.hartwig.hmftools.patientreporter.SampleMetadata;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.SampleReportFactory;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.pipeline.PipelineVersion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class AnalysedPatientReporter {

    private static final Logger LOGGER = LogManager.getLogger(AnalysedPatientReporter.class);

    @NotNull
    private final AnalysedReportData reportData;
    @NotNull
    private final String reportDate;

    public AnalysedPatientReporter(@NotNull final AnalysedReportData reportData, @NotNull final String reportDate) {
        this.reportData = reportData;
        this.reportDate = reportDate;
    }

    @NotNull
    public AnalysedPatientReport run(@NotNull SampleMetadata sampleMetadata, @NotNull PatientReporterConfig config) throws IOException {
        String patientId = reportData.limsModel().patientId(sampleMetadata.tumorSampleBarcode());
        PatientPrimaryTumor patientPrimaryTumor =
                PatientPrimaryTumorFunctions.findPrimaryTumorForPatient(reportData.patientPrimaryTumors(), patientId);

        SampleReport sampleReport = SampleReportFactory.fromLimsModel(sampleMetadata,
                reportData.limsModel(),
                patientPrimaryTumor,
                config.allowDefaultCohortConfig());

        String roseTcvFile = config.roseTsv();
        String clinicalSummary = config.addRose() && roseTcvFile != null ? RoseConclusionFile.read(roseTcvFile) : Strings.EMPTY;

        String specialRemark = reportData.specialRemarkModel().findSpecialRemarkForSample(sampleMetadata.tumorSampleId());

        String pipelineVersion = null;
        if (config.requirePipelineVersionFile()) {
            String pipelineVersionFile = config.pipelineVersionFile();
            assert pipelineVersionFile != null;
            pipelineVersion = PipelineVersionFile.majorDotMinorVersion(pipelineVersionFile);
            PipelineVersion.checkPipelineVersion(pipelineVersion, config.expectedPipelineVersion(), config.overridePipelineVersion());
        }

        GenomicAnalyzer genomicAnalyzer = new GenomicAnalyzer(reportData.germlineReportingModel());
        GenomicAnalysis genomicAnalysis = genomicAnalyzer.run(sampleMetadata.tumorSampleId(),
                sampleMetadata.refSampleId(),
                config,
                sampleReport.germlineReportingLevel(),
                reportData.knownFusionCache());

        GenomicAnalysis filteredAnalysis =
                ConsentFilterFunctions.filter(genomicAnalysis, sampleReport.germlineReportingLevel(), sampleReport.reportViralPresence());

        String qcForm = determineForNumber(genomicAnalysis.hasReliablePurity(), genomicAnalysis.impliedPurity());

        GenomicAnalysis overruledAnalysis = QualityOverruleFunctions.overrule(filteredAnalysis);
        GenomicAnalysis curateGeneName = CurationFunction.curation(overruledAnalysis);

        LOGGER.info("Loading CUPPA from {}", new File(config.cuppaResultCsv()).getParent());
        List<CuppaDataFile> cuppaEntries = CuppaDataFile.read(config.cuppaResultCsv());
        LOGGER.info(" Loaded {} entries from {}", cuppaEntries.size(), config.cuppaResultCsv());

        List<CuppaPrediction> predictions = CuppaPredictionFactory.create(cuppaEntries);
        CuppaPrediction best = predictions.get(0);
        CuppaReporting cuppaReporting = CuppaReportingFactory.createCuppaReportingData(best);

        LOGGER.info(" Predicted cancer type '{}' with likelihood {}", best.cancerType(), best.likelihood());

        List<PeachGenotype> peachGenotypes = curateGeneName.purpleQCStatus().contains(PurpleQCStatus.FAIL_CONTAMINATION)
                ? Lists.newArrayList()
                : loadPeachData(config.peachGenotypeTsv());

        List<PeachGenotype> peachGenotypesOverrule = sampleReport.reportPharmogenetics() ? peachGenotypes : Lists.newArrayList();

        Map<String, List<PeachGenotype>> peachMap = Maps.newHashMap();
        for (PeachGenotype peach : peachGenotypesOverrule) {
            if (peachMap.containsKey(peach.gene())) {
                List<PeachGenotype> curent = peachMap.get(peach.gene());
                curent.add(peach);
                peachMap.put(peach.gene(), curent);
            } else {
                peachMap.put(peach.gene(), Lists.newArrayList(peach));
            }
        }

        AnalysedPatientReport report = ImmutableAnalysedPatientReport.builder()
                .sampleReport(sampleReport)
                .qsFormNumber(qcForm)
                .clinicalSummary(clinicalSummary)
                .specialRemark(specialRemark)
                .pipelineVersion(pipelineVersion)
                .genomicAnalysis(curateGeneName)
                .cuppaReporting(
                        curateGeneName.purpleQCStatus().contains(PurpleQCStatus.FAIL_CONTAMINATION) || !curateGeneName.hasReliablePurity()
                                ? null
                                : cuppaReporting)
                .cuppaPlot(config.cuppaPlot())
                .circosPath(config.purpleCircosPlot())
                .comments(Optional.ofNullable(config.comments()))
                .isCorrectedReport(config.isCorrectedReport())
                .isCorrectedReportExtern(config.isCorrectedReportExtern())
                .signaturePath(reportData.signaturePath())
                .logoRVAPath(reportData.logoRVAPath())
                .logoCompanyPath(reportData.logoCompanyPath())
                .udiDi(reportData.udiDi())
                .peachGenotypes(peachMap)
                .reportDate(reportDate)
                .isWGSreport(true)
                .build();

        printReportState(report);

        return report;
    }

    @NotNull
    private static List<PeachGenotype> loadPeachData(@NotNull String peachGenotypeTsv) throws IOException {
        LOGGER.info("Loading peach genotypes from {}", new File(peachGenotypeTsv).getParent());
        List<PeachGenotype> peachGenotypes = PeachGenotypeFile.read(peachGenotypeTsv);
        LOGGER.info(" Loaded {} reportable genotypes from {}", peachGenotypes.size(), peachGenotypeTsv);
        return peachGenotypes;
    }

    @NotNull
    @VisibleForTesting
    static String determineForNumber(boolean hasReliablePurity, double purity) {
        return hasReliablePurity && purity > ReportResources.PURITY_CUTOFF
                ? QsFormNumber.FOR_080.display()
                : QsFormNumber.FOR_209.display();
    }

    private static void printReportState(@NotNull AnalysedPatientReport report) {
        LocalDate tumorArrivalDate = report.sampleReport().tumorArrivalDate();
        String formattedTumorArrivalDate =
                tumorArrivalDate != null ? DateTimeFormatter.ofPattern("dd-MMM-yyyy").format(tumorArrivalDate) : "N/A";

        LOGGER.info("Printing clinical and laboratory data for {}", report.sampleReport().tumorSampleId());
        LOGGER.info(" Tumor sample arrived at HMF on {}", formattedTumorArrivalDate);
        LOGGER.info(" Primary tumor details: {}{}",
                report.sampleReport().primaryTumorLocationString(),
                !report.sampleReport().primaryTumorTypeString().isEmpty()
                        ? " (" + report.sampleReport().primaryTumorTypeString() + ")"
                        : Strings.EMPTY);
        LOGGER.info(" Shallow seq purity: {}", report.sampleReport().shallowSeqPurityString());
        LOGGER.info(" Lab SOPs used: {}", report.sampleReport().labProcedures());
        LOGGER.info(" Clinical summary present: {}", (report.clinicalSummary() != null ? "yes" : "no"));
        LOGGER.info(" Special remark present: {}", (!report.specialRemark().isEmpty() ? "yes" : "no"));

        LOGGER.info(" Cohort: {}", report.sampleReport().cohort().cohortId());
        LOGGER.info(" Germline reporting level: {}", report.sampleReport().germlineReportingLevel());

        GenomicAnalysis analysis = report.genomicAnalysis();

        LOGGER.info("Printing genomic analysis results for {}:", report.sampleReport().tumorSampleId());
        if (report.cuppaReporting() != null) {
            LOGGER.info(" Molecular tissue origin conclusion: {}", report.cuppaReporting().interpretCancerType());
        }
        LOGGER.info(" Somatic variants to report: {}", analysis.reportableVariants().size());
        if (report.sampleReport().germlineReportingLevel() != LimsGermlineReportingLevel.NO_REPORTING) {
            LOGGER.info("  Number of variants known to exist in germline: {}", germlineOnly(analysis.reportableVariants()).size());
        } else {
            LOGGER.info("  Germline variants and evidence have been removed since no consent has been given");
        }
        LOGGER.info(" Number of gains and losses to report: {}", analysis.gainsAndLosses().size());
        LOGGER.info(" Gene fusions to report: {}", analysis.geneFusions().size());
        LOGGER.info(" Homozygous disruptions to report: {}", analysis.homozygousDisruptions().size());
        LOGGER.info(" Gene disruptions to report: {}", analysis.geneDisruptions().size());
        LOGGER.info(" Viruses to report: {}", analysis.reportableViruses().size());
        LOGGER.info(" Pharmacogenetics to report: {}", report.peachGenotypes().size());

        LOGGER.info(" CHORD analysis HRD prediction: {} ({})", analysis.chordHrdValue(), analysis.chordHrdStatus());
        LOGGER.info(" Microsatellite indels per Mb: {} ({})", analysis.microsatelliteIndelsPerMb(), analysis.microsatelliteStatus());
        LOGGER.info(" Tumor mutational load: {} ({})", analysis.tumorMutationalLoad(), analysis.tumorMutationalLoadStatus());
        LOGGER.info(" Tumor mutational burden: {}", analysis.tumorMutationalBurden());

        LOGGER.info("Printing actionability results for {}", report.sampleReport().tumorSampleId());
        LOGGER.info(" Tumor-specific evidence items found: {}", analysis.tumorSpecificEvidence().size());
        LOGGER.info(" Clinical trials matched to molecular profile: {}", analysis.clinicalTrials().size());
        LOGGER.info(" Off-label evidence items found: {}", analysis.offLabelEvidence().size());
    }

    @NotNull
    private static List<ReportableVariant> germlineOnly(@NotNull List<ReportableVariant> variants) {
        List<ReportableVariant> germlineOnly = Lists.newArrayList();
        for (ReportableVariant variant : variants) {
            if (variant.source() == ReportableVariantSource.GERMLINE) {
                germlineOnly.add(variant);
            }
        }
        return germlineOnly;
    }
}