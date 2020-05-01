package com.hartwig.hmftools.patientreporter;

import java.io.IOException;
import java.text.DecimalFormat;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ChordFileReader;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocationFunctions;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingChoice;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.qc.PurpleQC;
import com.hartwig.hmftools.common.purple.qc.PurpleQCFile;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableDisruption;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableDisruptionFile;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusionFile;
import com.hartwig.hmftools.patientreporter.actionability.ClinicalTrialFactory;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberAnalysis;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberAnalyzer;
import com.hartwig.hmftools.patientreporter.homozygousdisruption.HomozygousDisruptionAnalyzer;
import com.hartwig.hmftools.patientreporter.homozygousdisruption.ReportableHomozygousDisruption;
import com.hartwig.hmftools.patientreporter.structural.SvAnalysis;
import com.hartwig.hmftools.patientreporter.structural.SvAnalyzer;
import com.hartwig.hmftools.patientreporter.variants.ReportVariantAnalysis;
import com.hartwig.hmftools.patientreporter.variants.ReportableGermlineVariant;
import com.hartwig.hmftools.patientreporter.variants.ReportableVariant;
import com.hartwig.hmftools.patientreporter.variants.ReportableVariantAnalyzer;
import com.hartwig.hmftools.patientreporter.variants.germline.BachelorFile;
import com.hartwig.hmftools.patientreporter.variants.germline.FilterGermlineVariants;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineVariant;
import com.hartwig.hmftools.patientreporter.variants.somatic.SomaticVariantAnalysis;
import com.hartwig.hmftools.patientreporter.variants.somatic.SomaticVariantAnalyzer;
import com.hartwig.hmftools.patientreporter.viralInsertion.ViralInsertion;
import com.hartwig.hmftools.patientreporter.viralInsertion.ViralInsertionAnalyzer;

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
            @NotNull String purpleGeneCnvTsv, @NotNull String somaticVariantVcf, @NotNull String bachelorTsv, @NotNull String linxFusionTsv,
            @NotNull String linxDisruptionTsv, @NotNull String linxViralInsertionTsv, @NotNull String linxDriversTsv,
            @NotNull String chordPredictionTxt, @NotNull String circosFile, @Nullable String comments, boolean correctedReport,
            boolean unofficialReport) throws IOException {
        PatientTumorLocation patientTumorLocation =
                PatientTumorLocationFunctions.findPatientTumorLocationForSample(reportData.patientTumorLocations(),
                        sampleMetadata.tumorSampleId());

        SampleReport sampleReport = SampleReportFactory.fromLimsAndHospitalModel(sampleMetadata,
                reportData.limsModel(),
                reportData.limsWideModel(),
                reportData.hospitalModel(),
                patientTumorLocation);

        CopyNumberAnalysis copyNumberAnalysis = analyzeCopyNumbers(purplePurityTsv, purpleQCFile, purpleGeneCnvTsv, patientTumorLocation);
        SomaticVariantAnalysis somaticVariantAnalysis =
                analyzeSomaticVariants(sampleMetadata.tumorSampleId(), somaticVariantVcf, copyNumberAnalysis.exomeGeneCopyNumbers());

        ChordAnalysis chordAnalysis = analyzeChord(chordPredictionTxt);
        List<ReportableGermlineVariant> germlineVariantsToReport = analyzeGermlineVariants(sampleMetadata.tumorSampleBarcode(),
                bachelorTsv,
                copyNumberAnalysis,
                somaticVariantAnalysis,
                chordAnalysis);

        ReportVariantAnalysis reportableVariantsAnalysis =
                ReportableVariantAnalyzer.mergeSomaticAndGermlineVariants(somaticVariantAnalysis.variantsToReport(),
                        somaticVariantAnalysis.driverCatalog(),
                        reportData.driverGeneView(),
                        germlineVariantsToReport,
                        reportData.germlineReportingModel(),
                        reportData.limsModel().germlineReportingChoice(sampleMetadata.tumorSampleBarcode()),
                        reportData.actionabilityAnalyzer(),
                        patientTumorLocation);

        SvAnalysis svAnalysis = analyzeStructuralVariants(linxFusionTsv, linxDisruptionTsv, patientTumorLocation);
        List<ReportableHomozygousDisruption> reportableHomozygousDisruptions = extractHomozygousDisruptionsFromLinxDrivers(linxDriversTsv);
        List<ViralInsertion> viralInsertions = analyzeViralInsertions(linxViralInsertionTsv,
                reportData.limsModel().viralInsertionChoice(sampleMetadata.tumorSampleBarcode(), sampleMetadata.tumorSampleId()));

        String clinicalSummary = reportData.summaryModel().findSummaryForSample(sampleMetadata.tumorSampleId());

        List<EvidenceItem> allEvidenceItems = Lists.newArrayList();
        allEvidenceItems.addAll(reportableVariantsAnalysis.evidenceItems());
        allEvidenceItems.addAll(copyNumberAnalysis.evidenceItems());
        allEvidenceItems.addAll(svAnalysis.evidenceItems());

        List<EvidenceItem> nonTrials = ReportableEvidenceItemFactory.extractNonTrials(allEvidenceItems);
        AnalysedPatientReport report = ImmutableAnalysedPatientReport.builder()
                .sampleReport(sampleReport)
                .impliedPurity(copyNumberAnalysis.purity())
                .hasReliablePurity(copyNumberAnalysis.hasReliablePurity())
                .hasReliableQuality(copyNumberAnalysis.hasReliableQuality())
                .averageTumorPloidy(copyNumberAnalysis.ploidy())
                .clinicalSummary(clinicalSummary)
                .tumorSpecificEvidence(nonTrials.stream().filter(EvidenceItem::isOnLabel).collect(Collectors.toList()))
                .clinicalTrials(ClinicalTrialFactory.extractOnLabelTrials(allEvidenceItems))
                .offLabelEvidence(nonTrials.stream().filter(item -> !item.isOnLabel()).collect(Collectors.toList()))
                .reportableVariants(reportableVariantsAnalysis.variantsToReport())
                .microsatelliteIndelsPerMb(copyNumberAnalysis.microsatelliteIndelsPerMb())
                .microsatelliteStatus(copyNumberAnalysis.microsatelliteStatus())
                .tumorMutationalLoad((int) Math.round(copyNumberAnalysis.tumorMutationalLoad()))
                .tumorMutationalLoadStatus(copyNumberAnalysis.tumorMutationalLoadStatus())
                .tumorMutationalBurden(copyNumberAnalysis.tumorMutationalBurdenPerMb())
                .chordAnalysis(chordAnalysis)
                .hrdStatus(ChordStatus.fromHRD(chordAnalysis.hrdValue()))
                .gainsAndLosses(copyNumberAnalysis.reportableGainsAndLosses())
                .geneFusions(svAnalysis.reportableFusions())
                .geneDisruptions(svAnalysis.reportableDisruptions())
                .reportableHomozygousDisruptions(reportableHomozygousDisruptions)
                .viralInsertions(viralInsertions)
                .circosPath(circosFile)
                .comments(Optional.ofNullable(comments))
                .isCorrectedReport(correctedReport)
                .isUnofficialReport(unofficialReport)
                .signaturePath(reportData.signaturePath())
                .logoRVAPath(reportData.logoRVAPath())
                .logoCompanyPath(reportData.logoCompanyPath())
                .build();

        printReportState(report);

        return report;
    }

    @NotNull
    private static List<ReportableHomozygousDisruption> extractHomozygousDisruptionsFromLinxDrivers(@NotNull String linxDriversTsv)
            throws IOException {
        return HomozygousDisruptionAnalyzer.extractFromLinxDriversTsv(linxDriversTsv);
    }

    @Nullable
    private static List<ViralInsertion> analyzeViralInsertions(@NotNull String linxViralInsertionTsv, boolean viralInsertionsChoice)
            throws IOException {
        return ViralInsertionAnalyzer.loadViralInsertions(linxViralInsertionTsv, viralInsertionsChoice);
    }

    @NotNull
    private CopyNumberAnalysis analyzeCopyNumbers(@NotNull String purplePurityTsv, @NotNull String purpleQCFile,
            @NotNull String purpleGeneCnvTsv, @Nullable PatientTumorLocation patientTumorLocation) throws IOException {
        PurityContext purityContext = FittedPurityFile.read(purplePurityTsv);
        PurpleQC purpleQC = PurpleQCFile.read(purpleQCFile);

        LOGGER.info("Loaded purple sample data from {}", purplePurityTsv);
        LOGGER.info(" Purple purity: {}", new DecimalFormat("#'%'").format(purityContext.bestFit().purity() * 100));
        LOGGER.info(" Purple average tumor ploidy: {}", purityContext.bestFit().ploidy());
        LOGGER.info(" Purple status: {}", purityContext.status());
        LOGGER.info(" WGD happened: {}", purityContext.wholeGenomeDuplication() ? "yes" : "no");
        LOGGER.info("Loaded purple QC data from {}", purpleQCFile);
        LOGGER.info(" Purple QC status: {}", purpleQC.status());

        List<GeneCopyNumber> exomeGeneCopyNumbers = GeneCopyNumberFile.read(purpleGeneCnvTsv);
        LOGGER.info("Loaded {} gene copy numbers from {}", exomeGeneCopyNumbers.size(), purpleGeneCnvTsv);

        return CopyNumberAnalyzer.run(purityContext,
                purpleQC,
                exomeGeneCopyNumbers,
                reportData.actionabilityAnalyzer(),
                patientTumorLocation);
    }

    @NotNull
    private SomaticVariantAnalysis analyzeSomaticVariants(@NotNull String sample, @NotNull String somaticVariantVcf,
            @NotNull List<GeneCopyNumber> exomeGeneCopyNumbers) throws IOException {
        List<SomaticVariant> variants = SomaticVariantFactory.passOnlyInstance().fromVCFFile(sample, somaticVariantVcf);
        LOGGER.info("Loaded {} PASS somatic variants from {}", variants.size(), somaticVariantVcf);

        return SomaticVariantAnalyzer.run(variants, reportData.driverGeneView(), exomeGeneCopyNumbers);
    }

    @NotNull
    private List<ReportableGermlineVariant> analyzeGermlineVariants(@NotNull String sampleBarcode, @NotNull String bachelorTsv,
            @NotNull CopyNumberAnalysis copyNumberAnalysis, @NotNull SomaticVariantAnalysis somaticVariantAnalysis,
            @NotNull ChordAnalysis chordAnalysis) throws IOException {
        List<GermlineVariant> variants =
                BachelorFile.loadBachelorTsv(bachelorTsv).stream().filter(GermlineVariant::passFilter).collect(Collectors.toList());
        LOGGER.info("Loaded {} PASS germline variants from {}", variants.size(), bachelorTsv);

        LimsGermlineReportingChoice germlineChoice = reportData.limsModel().germlineReportingChoice(sampleBarcode);
        if (germlineChoice != LimsGermlineReportingChoice.NO_REPORTING) {
            LOGGER.info(" Patient has given the following germline consent: '{}'", germlineChoice);
            return FilterGermlineVariants.filterGermlineVariantsForReporting(variants,
                    reportData.driverGeneView(),
                    reportData.germlineReportingModel(),
                    copyNumberAnalysis.exomeGeneCopyNumbers(),
                    somaticVariantAnalysis.variantsToReport(),
                    chordAnalysis);
        } else {
            LOGGER.info(" No consent has been given for germline reporting. No germline variants will be reported!");
            return Lists.newArrayList();
        }
    }

    @NotNull
    private SvAnalysis analyzeStructuralVariants(@NotNull String linxFusionTsv, @NotNull String linxDisruptionTsv,
            @Nullable PatientTumorLocation patientTumorLocation) throws IOException {
        List<ReportableGeneFusion> fusions = ReportableGeneFusionFile.read(linxFusionTsv);
        LOGGER.info("Loaded {} fusions from {}", fusions.size(), linxFusionTsv);

        List<ReportableDisruption> disruptions = ReportableDisruptionFile.read(linxDisruptionTsv);
        LOGGER.info("Loaded {} disruptions from {}", disruptions.size(), linxDisruptionTsv);

        return SvAnalyzer.run(fusions, disruptions, reportData.actionabilityAnalyzer(), patientTumorLocation);
    }

    @NotNull
    private static ChordAnalysis analyzeChord(@NotNull String chordPredictionTxt) throws IOException {
        ChordAnalysis chord = ChordFileReader.read(chordPredictionTxt);
        LOGGER.info("Loaded CHORD analysis from {}", chordPredictionTxt);
        return chord;
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
        LOGGER.info(" Shallow seq purity: {}", report.sampleReport().purityShallowSeq());
        LOGGER.info(" Lab SOPs used: {}", report.sampleReport().labProcedures());
        LOGGER.info(" Clinical summary present: {}", (!report.clinicalSummary().isEmpty() ? "yes" : "no"));

        List<ReportableVariant> variantsWithNotify =
                report.reportableVariants().stream().filter(ReportableVariant::notifyClinicalGeneticist).collect(Collectors.toList());
        LOGGER.info("Printing genomic analysis results for {}:", report.sampleReport().tumorSampleId());
        LOGGER.info(" Somatic variants to report: {}", report.reportableVariants().size());
        LOGGER.info("  Variants for which to notify clinical geneticist: {}", variantsWithNotify.size());
        LOGGER.info(" Microsatellite Indels per Mb: {}", report.microsatelliteIndelsPerMb());
        LOGGER.info(" Tumor mutational load: {}", report.tumorMutationalLoad());
        LOGGER.info(" Tumor mutational burden: {}", report.tumorMutationalBurden());
        LOGGER.info(" CHORD analysis HRD prediction: {}", report.chordAnalysis().hrdValue());
        LOGGER.info(" Number of gains and losses to report: {}", report.gainsAndLosses().size());
        LOGGER.info(" Gene fusions to report : {}", report.geneFusions().size());
        LOGGER.info(" Gene disruptions to report : {}", report.geneDisruptions().size());

        LOGGER.info("Printing actionability results for {}", report.sampleReport().tumorSampleId());
        LOGGER.info(" Tumor-specific evidence items found: {}", report.tumorSpecificEvidence().size());
        LOGGER.info(" Off-label evidence items found: {}", report.offLabelEvidence().size());
        LOGGER.info(" Clinical trials matched to molecular profile: {}", report.clinicalTrials().size());
    }
}
