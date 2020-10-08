package com.hartwig.hmftools.patientreporter;

import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.clinical.PatientTumorLocation;
import com.hartwig.hmftools.common.clinical.PatientTumorLocationFunctions;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.protect.GenomicAnalysis;
import com.hartwig.hmftools.protect.GenomicAnalyzer;
import com.hartwig.hmftools.protect.variants.ReportableVariant;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class AnalysedPatientReporter {

    private static final Logger LOGGER = LogManager.getLogger(AnalysedPatientReporter.class);

    @NotNull
    private final AnalysedReportData reportData;

    AnalysedPatientReporter(@NotNull final AnalysedReportData reportData) {
        this.reportData = reportData;
    }

    @NotNull
    AnalysedPatientReport run(@NotNull SampleMetadata sampleMetadata, @NotNull String purplePurityTsv, @NotNull String purpleQCFile,
            @NotNull String purpleGeneCnvTsv, @NotNull String purpleDriverCatalogTsv, @NotNull String somaticVariantVcf,
            @NotNull String bachelorTsv, @NotNull String linxFusionTsv, @NotNull String linxDisruptionTsv,
            @NotNull String linxViralInsertionTsv, @NotNull String linxDriversTsv, @NotNull String chordPredictionTxt,
            @NotNull String circosFile, @Nullable String comments, boolean correctedReport) throws IOException {
        PatientTumorLocation patientTumorLocation =
                PatientTumorLocationFunctions.findTumorLocationForSample(reportData.patientTumorLocations(),
                        sampleMetadata.tumorSampleId());

        SampleReport sampleReport = SampleReportFactory.fromLimsModel(sampleMetadata, reportData.limsModel(), patientTumorLocation);
        LimsGermlineReportingLevel germlineChoice = reportData.limsModel().germlineReportingChoice(sampleMetadata.tumorSampleBarcode());
        boolean reportViralInsertions = reportData.limsModel().reportViralInsertions(sampleMetadata.tumorSampleBarcode());

        GenomicAnalyzer genomicAnalyzer = new GenomicAnalyzer(reportData.actionabilityAnalyzer(), reportData.germlineReportingModel());
        GenomicAnalysis genomicAnalysis = genomicAnalyzer.analyze(sampleMetadata.tumorSampleId(),
                patientTumorLocation,
                germlineChoice,
                reportViralInsertions,
                purplePurityTsv,
                purpleQCFile,
                purpleGeneCnvTsv,
                purpleDriverCatalogTsv,
                somaticVariantVcf,
                bachelorTsv,
                linxFusionTsv,
                linxDisruptionTsv,
                linxViralInsertionTsv,
                linxDriversTsv,
                chordPredictionTxt);

        String clinicalSummary = reportData.summaryModel().findSummaryForSample(sampleMetadata.tumorSampleId());

        AnalysedPatientReport report = ImmutableAnalysedPatientReport.builder()
                .sampleReport(sampleReport)
                .qsFormNumber(determineForNumber(genomicAnalysis.hasReliablePurity(), genomicAnalysis.impliedPurity()))
                .clinicalSummary(clinicalSummary)
                .genomicAnalysis(genomicAnalysis)
                .circosPath(circosFile)
                .comments(Optional.ofNullable(comments))
                .isCorrectedReport(correctedReport)
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
        LOGGER.info(" Primary tumor location: {}{}",
                report.sampleReport().primaryTumorLocationString(),
                !report.sampleReport().cancerSubTypeString().isEmpty()
                        ? " (" + report.sampleReport().cancerSubTypeString() + ")"
                        : Strings.EMPTY);
        LOGGER.info(" Shallow seq purity: {}", report.sampleReport().shallowSeqPurityString());
        LOGGER.info(" Lab SOPs used: {}", report.sampleReport().labProcedures());
        LOGGER.info(" Clinical summary present: {}", (!report.clinicalSummary().isEmpty() ? "yes" : "no"));

        GenomicAnalysis analysis = report.genomicAnalysis();
        List<ReportableVariant> variantsWithNotify =
                analysis.reportableVariants().stream().filter(ReportableVariant::notifyClinicalGeneticist).collect(Collectors.toList());
        LOGGER.info("Printing genomic analysis results for {}:", report.sampleReport().tumorSampleId());
        LOGGER.info(" Somatic variants to report: {}", analysis.reportableVariants().size());
        LOGGER.info("  Variants for which to notify clinical geneticist: {}", variantsWithNotify.size());
        LOGGER.info(" Microsatellite indels per Mb: {} ({})", analysis.microsatelliteIndelsPerMb(), analysis.microsatelliteStatus());
        LOGGER.info(" Tumor mutational load: {} ({})", analysis.tumorMutationalLoad(), analysis.tumorMutationalLoadStatus());
        LOGGER.info(" Tumor mutational burden: {}", analysis.tumorMutationalBurden());
        LOGGER.info(" CHORD analysis HRD prediction: {} ({})", analysis.chordHrdValue(), analysis.chordHrdStatus());
        LOGGER.info(" Number of gains and losses to report: {}", analysis.gainsAndLosses().size());
        LOGGER.info(" Gene fusions to report: {}", analysis.geneFusions().size());
        LOGGER.info(" Gene disruptions to report: {}", analysis.geneDisruptions().size());
        LOGGER.info(" Viral insertions to report: {}", analysis.viralInsertions() != null ? analysis.viralInsertions().size() : "0");

        LOGGER.info("Printing actionability results for {}", report.sampleReport().tumorSampleId());
        LOGGER.info(" Tumor-specific evidence items found: {}", analysis.tumorSpecificEvidence().size());
        LOGGER.info(" Clinical trials matched to molecular profile: {}", analysis.clinicalTrials().size());
        LOGGER.info(" Off-label evidence items found: {}", analysis.offLabelEvidence().size());
    }
}
