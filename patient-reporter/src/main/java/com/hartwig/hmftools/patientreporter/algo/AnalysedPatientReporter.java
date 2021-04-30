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
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.runcontext.MetaDataResolver;
import com.hartwig.hmftools.patientreporter.QsFormNumber;
import com.hartwig.hmftools.patientreporter.SampleMetadata;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.SampleReportFactory;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cuppa.ImmutableMolecularTissueOrigin;
import com.hartwig.hmftools.patientreporter.cuppa.MolecularTissueOrigin;
import com.hartwig.hmftools.patientreporter.cuppa.MolecularTissueOriginFactory;
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
    public AnalysedPatientReport run(@NotNull SampleMetadata sampleMetadata, @NotNull String purplePurityTsv, @NotNull String purpleQCFile,
            @NotNull String purpleDriverCatalogSomaticTsv, @NotNull String purpleDriverCatalogGermlineTsv,
            @NotNull String purpleSomaticVariantVcf, @NotNull String purpleGermlineVariantVcf, @NotNull String linxFusionTsv,
            @NotNull String linxBreakendTsv, @NotNull String linxDriversTsv,
            @NotNull String chordPredictionTxt, @NotNull String circosFile, @NotNull String protectEvidenceTsv, @Nullable String comments,
            boolean correctedReport, @NotNull String pipelineVersionFile, @NotNull String molecularTissueOriginTsv,
            @NotNull String molecularTissueOriginPlot, @NotNull String virusBreakendTsv) throws IOException {
        // TODO Specific COLO handling doesn't belong in patient reporter!
        String patientId = sampleMetadata.patientId().startsWith("COLO829") ? "COLO829" : sampleMetadata.patientId();
        PatientPrimaryTumor patientPrimaryTumor =
                PatientPrimaryTumorFunctions.findPrimaryTumorForPatient(reportData.patientPrimaryTumors(), patientId);

        SampleReport sampleReport = SampleReportFactory.fromLimsModel(sampleMetadata, reportData.limsModel(), patientPrimaryTumor);

        GenomicAnalyzer genomicAnalyzer = new GenomicAnalyzer();
        GenomicAnalysis genomicAnalysis = genomicAnalyzer.run(sampleMetadata.tumorSampleId(),
                purplePurityTsv,
                purpleQCFile,
                purpleDriverCatalogSomaticTsv,
                purpleDriverCatalogGermlineTsv,
                purpleSomaticVariantVcf,
                purpleGermlineVariantVcf,
                linxFusionTsv,
                linxBreakendTsv,
                linxDriversTsv,
                chordPredictionTxt,
                protectEvidenceTsv, virusBreakendTsv);

        ConsentFilterFunctions consentFilterFunctions = new ConsentFilterFunctions();

        GenomicAnalysis filteredAnalysis = consentFilterFunctions.filterAndOverruleForConsent(genomicAnalysis,
                sampleReport.germlineReportingLevel(),
                sampleReport.reportViralInsertions());

        String clinicalSummary = reportData.summaryModel().findSummaryForSample(sampleMetadata.tumorSampleId(), sampleReport.cohort());

        // TODO Move format to PDF. Make pipeline version in report @Nullable.
        String pipelineVersion = !pipelineVersionFile.isEmpty()
                ? MetaDataResolver.majorDotMinorVersion(new File(pipelineVersionFile))
                : "No pipeline version is known";

        MolecularTissueOrigin molecularTissueOrigin = ImmutableMolecularTissueOrigin.builder()
                .molecularTissueOriginResult(MolecularTissueOriginFactory.readMolecularTissueOriginResult(molecularTissueOriginTsv))
                .molecularTissueOriginPlot(molecularTissueOriginPlot)
                .build();

        AnalysedPatientReport report = ImmutableAnalysedPatientReport.builder()
                .sampleReport(sampleReport)
                .qsFormNumber(determineForNumber(genomicAnalysis.hasReliablePurity(), genomicAnalysis.impliedPurity()))
                .clinicalSummary(clinicalSummary)
                .pipelineVersion(pipelineVersion)
                .genomicAnalysis(filteredAnalysis)
                .molecularTissueOrigin(molecularTissueOrigin)
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
        LOGGER.info(" Primary tumor details: {}{}",
                report.sampleReport().primaryTumorLocationString(),
                !report.sampleReport().primaryTumorTypeString().isEmpty()
                        ? " (" + report.sampleReport().primaryTumorTypeString() + ")"
                        : Strings.EMPTY);
        LOGGER.info(" Shallow seq purity: {}", report.sampleReport().shallowSeqPurityString());
        LOGGER.info(" Lab SOPs used: {}", report.sampleReport().labProcedures());
        LOGGER.info(" Clinical summary present: {}", (!report.clinicalSummary().isEmpty() ? "yes" : "no"));

        GenomicAnalysis analysis = report.genomicAnalysis();

        LOGGER.info("Printing genomic analysis results for {}:", report.sampleReport().tumorSampleId());
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
