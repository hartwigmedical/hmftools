package com.hartwig.hmftools.patientreporter;

import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;
import java.util.Optional;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.ClinicalTrial;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.clinical.PatientTumorLocation;
import com.hartwig.hmftools.common.clinical.PatientTumorLocationFunctions;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.data.EventFilter;
import com.hartwig.hmftools.protect.GenomicAnalysis;
import com.hartwig.hmftools.protect.GenomicAnalyzer;
import com.hartwig.hmftools.protect.variants.ReportableVariant;
import com.hartwig.hmftools.protect.variants.ReportableVariantSource;

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
            @NotNull String purpleDriverCatalogTsv, @NotNull String purpleSomaticVariantVcf, @NotNull String bachelorTsv,
            @NotNull String linxFusionTsv, @NotNull String linxBreakendTsv, @NotNull String linxViralInsertionTsv,
            @NotNull String linxDriversTsv, @NotNull String chordPredictionTxt, @NotNull String circosFile, @Nullable String comments,
            boolean correctedReport) throws IOException {
        PatientTumorLocation patientTumorLocation =
                PatientTumorLocationFunctions.findTumorLocationForSample(reportData.patientTumorLocations(),
                        sampleMetadata.tumorSampleId());

        SampleReport sampleReport = SampleReportFactory.fromLimsModel(sampleMetadata, reportData.limsModel(), patientTumorLocation);

        GenomicAnalyzer genomicAnalyzer = new GenomicAnalyzer(reportData.actionabilityAnalyzer(), reportData.germlineReportingModel());
        GenomicAnalysis genomicAnalysis = genomicAnalyzer.run(sampleMetadata.tumorSampleId(),
                patientTumorLocation,
                purplePurityTsv,
                purpleQCFile,
                purpleDriverCatalogTsv,
                purpleSomaticVariantVcf,
                bachelorTsv,
                linxFusionTsv,
                linxBreakendTsv,
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
                !report.sampleReport().primaryTumorTypeString().isEmpty()
                        ? " (" + report.sampleReport().primaryTumorTypeString() + ")"
                        : Strings.EMPTY);
        LOGGER.info(" Shallow seq purity: {}", report.sampleReport().shallowSeqPurityString());
        LOGGER.info(" Lab SOPs used: {}", report.sampleReport().labProcedures());
        LOGGER.info(" Clinical summary present: {}", (!report.clinicalSummary().isEmpty() ? "yes" : "no"));

        GenomicAnalysis analysis = report.genomicAnalysis();

        List<ClinicalTrial> filteredClinicalTrials = EventFilter.removeEvidenceOnFilteredGermlineVariants(analysis.clinicalTrials(),
                analysis.reportableVariants(),
                report.sampleReport().germlineReportingLevel());

        List<EvidenceItem> filteredTumorSpecificEvidence = EventFilter.removeEvidenceOnFilteredGermlineVariants(analysis.tumorSpecificEvidence(),
                analysis.reportableVariants(),
                report.sampleReport().germlineReportingLevel());

        List<EvidenceItem> filteredOffLabelEvidence = EventFilter.removeEvidenceOnFilteredGermlineVariants(analysis.offLabelEvidence(),
                analysis.reportableVariants(),
                report.sampleReport().germlineReportingLevel());

        LOGGER.info("Printing genomic analysis results for {}:", report.sampleReport().tumorSampleId());
        LOGGER.info(" Somatic variants to report: {}", analysis.reportableVariants().size());
        if (report.sampleReport().germlineReportingLevel() != LimsGermlineReportingLevel.NO_REPORTING) {
            LOGGER.info("  Number of variants also present in germline: {}", germlineOnly(analysis.reportableVariants()).size());
        } else {
            LOGGER.info("  Germline variants and evidence will be removed since no consent has been given");
        }
        LOGGER.info(" Number of gains and losses to report: {}", analysis.gainsAndLosses().size());
        LOGGER.info(" Gene fusions to report: {}", analysis.geneFusions().size());
        LOGGER.info(" Homozygous disruptions to report: {}", analysis.homozygousDisruptions().size());
        LOGGER.info(" Gene disruptions to report: {}", analysis.geneDisruptions().size());
        LOGGER.info(" Viral insertions to report: {}", analysis.viralInsertions().size());
        LOGGER.info(" CHORD analysis HRD prediction: {} ({})", analysis.chordHrdValue(), analysis.chordHrdStatus());
        LOGGER.info(" Microsatellite indels per Mb: {} ({})", analysis.microsatelliteIndelsPerMb(), analysis.microsatelliteStatus());
        LOGGER.info(" Tumor mutational load: {} ({})", analysis.tumorMutationalLoad(), analysis.tumorMutationalLoadStatus());
        LOGGER.info(" Tumor mutational burden: {}", analysis.tumorMutationalBurden());

        LOGGER.info("Printing actionability results for {}", report.sampleReport().tumorSampleId());
        LOGGER.info(" Tumor-specific evidence items found: {}", filteredTumorSpecificEvidence.size());
        LOGGER.info(" Clinical trials matched to molecular profile: {}", filteredClinicalTrials.size());
        LOGGER.info(" Off-label evidence items found: {}", filteredOffLabelEvidence.size());
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
