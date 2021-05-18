package com.hartwig.hmftools.patientreporter.algo;

import java.io.File;
import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;
import java.util.Optional;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.clinical.PatientPrimaryTumor;
import com.hartwig.hmftools.common.clinical.PatientPrimaryTumorFunctions;
import com.hartwig.hmftools.common.cuppa.ImmutableMolecularTissueOrigin;
import com.hartwig.hmftools.common.cuppa.MolecularTissueOrigin;
import com.hartwig.hmftools.common.cuppa.MolecularTissueOriginFile;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.runcontext.MetaDataResolver;
import com.hartwig.hmftools.patientreporter.PatientReporterConfig;
import com.hartwig.hmftools.patientreporter.QsFormNumber;
import com.hartwig.hmftools.patientreporter.SampleMetadata;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.SampleReportFactory;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.protect.purple.ReportableVariant;
import com.hartwig.hmftools.protect.purple.ReportableVariantSource;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class AnalysedPatientReporter {

    private static final Logger LOGGER = LogManager.getLogger(AnalysedPatientReporter.class);

    @NotNull
    private final AnalysedReportData reportData;

    public AnalysedPatientReporter(@NotNull final AnalysedReportData reportData) {
        this.reportData = reportData;
    }

    @NotNull
    public AnalysedPatientReport run(@NotNull SampleMetadata sampleMetadata, @NotNull PatientReporterConfig config) throws IOException {
        String patientId = reportData.limsModel().patientId(sampleMetadata.tumorSampleBarcode());
        PatientPrimaryTumor patientPrimaryTumor =
                PatientPrimaryTumorFunctions.findPrimaryTumorForPatient(reportData.patientPrimaryTumors(), patientId);

        SampleReport sampleReport = SampleReportFactory.fromLimsModel(sampleMetadata, reportData.limsModel(), patientPrimaryTumor);

        String clinicalSummary = reportData.summaryModel().findSummaryForSample(sampleMetadata.tumorSampleId(), sampleReport.cohort());

        String pipelineVersion = MetaDataResolver.majorDotMinorVersion(new File(config.pipelineVersionFile()));
        checkPipelineVersion(pipelineVersion, config);

        GenomicAnalyzer genomicAnalyzer = new GenomicAnalyzer(reportData.germlineReportingModel(),
                reportData.taxonomyDb(),
                reportData.virusInterpretationModel(),
                reportData.virusBlackListModel());
        GenomicAnalysis genomicAnalysis =
                genomicAnalyzer.run(sampleMetadata.tumorSampleId(), config, sampleReport.germlineReportingLevel());

        GenomicAnalysis filteredAnalysis = ConsentFilterFunctions.filter(genomicAnalysis,
                sampleReport.germlineReportingLevel(),
                sampleReport.reportViralInsertions(),
                sampleReport.cohort().reportPeach());

        GenomicAnalysis overruledAnalysis = QualityOverruleFunctions.overrule(filteredAnalysis);

        LOGGER.info("Loading CUPPA result from {}", new File(config.molecularTissueOriginTxt()).getParent());
        MolecularTissueOrigin molecularTissueOrigin = ImmutableMolecularTissueOrigin.builder()
                .conclusion(MolecularTissueOriginFile.read(config.molecularTissueOriginTxt()))
                .plotPath(config.molecularTissueOriginPlot())
                .build();
        LOGGER.info(" Molecular tissue origin conclusion: {}", molecularTissueOrigin.conclusion());

        AnalysedPatientReport report = ImmutableAnalysedPatientReport.builder()
                .sampleReport(sampleReport)
                .qsFormNumber(determineForNumber(genomicAnalysis.hasReliablePurity(), genomicAnalysis.impliedPurity()))
                .clinicalSummary(clinicalSummary)
                .pipelineVersion(pipelineVersion)
                .genomicAnalysis(overruledAnalysis)
                .molecularTissueOrigin(molecularTissueOrigin)
                .circosPath(config.purpleCircosPlot())
                .comments(Optional.ofNullable(config.comments()))
                .isCorrectedReport(config.isCorrectedReport())
                .signaturePath(reportData.signaturePath())
                .logoRVAPath(reportData.logoRVAPath())
                .logoCompanyPath(reportData.logoCompanyPath())
                .build();

        printReportState(report);

        return report;
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
        LOGGER.info(" Clinical summary present: {}", (!report.clinicalSummary().isEmpty() ? "yes" : "no"));
        LOGGER.info(" Cohort: {}", report.sampleReport().cohort().cohortId());
        LOGGER.info(" Germline reporting level: {}", report.sampleReport().germlineReportingLevel());

        GenomicAnalysis analysis = report.genomicAnalysis();

        LOGGER.info("Printing genomic analysis results for {}:", report.sampleReport().tumorSampleId());
        LOGGER.info(" Molecular tissue origin conclusion: {}", report.molecularTissueOrigin().conclusion());
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
        LOGGER.info(" Virus breakend to report: {}", analysis.virusBreakends().size());
        LOGGER.info(" Pharmacogenetics to report: {}", analysis.peachGenotypes().size());

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

    private static void checkPipelineVersion(@Nullable String pipelineVersion, @NotNull PatientReporterConfig config) {
        if (!config.overridePipelineVersion()) {
            if (pipelineVersion != null && !pipelineVersion.equals(config.expectedPipelineVersion())) {
                LOGGER.warn("The expected pipeline version {} is different than the real pipeline version {}!",
                        pipelineVersion,
                        config.expectedPipelineVersion());
            }
        } else if (config.overridePipelineVersion()) {
            LOGGER.warn("Pipeline version is overridden!");
        }
    }
}
