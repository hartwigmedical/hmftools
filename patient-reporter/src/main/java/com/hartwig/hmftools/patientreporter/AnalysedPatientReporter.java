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
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocationFunctions;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingChoice;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.enrich.VariantContextEnrichmentFactory;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableDisruption;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableDisruptionFile;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusionFile;
import com.hartwig.hmftools.patientreporter.actionability.ClinicalTrialFactory;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberAnalysis;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberAnalyzer;
import com.hartwig.hmftools.patientreporter.structural.SvAnalysis;
import com.hartwig.hmftools.patientreporter.structural.SvAnalyzer;
import com.hartwig.hmftools.patientreporter.variants.ReportableVariant;
import com.hartwig.hmftools.patientreporter.variants.ReportableVariantAnalyzer;
import com.hartwig.hmftools.patientreporter.variants.germline.BachelorFile;
import com.hartwig.hmftools.patientreporter.variants.germline.FilterGermlineVariants;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineVariant;
import com.hartwig.hmftools.patientreporter.variants.somatic.SomaticVariantAnalysis;
import com.hartwig.hmftools.patientreporter.variants.somatic.SomaticVariantAnalyzer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
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
    public AnalysedPatientReport run(@NotNull String tumorSample, @NotNull String refSample, @NotNull String purplePurityTsv,
            @NotNull String purpleGeneCnvTsv, @NotNull String somaticVariantVcf, @NotNull String linxFusionTsv,
            @NotNull String linxDisruptionTsv, @Nullable String bachelorCsv, @NotNull String chordPredictionFile,
            @NotNull String circosFile, @Nullable String comments, @Nullable String correctTitle) throws IOException {
        PatientTumorLocation patientTumorLocation =
                PatientTumorLocationFunctions.findPatientTumorLocationForSample(reportData.patientTumorLocations(), tumorSample);

        SampleReport sampleReport = SampleReportFactory.fromLimsAndHospitalModel(tumorSample,
                refSample,
                reportData.limsModel(),
                reportData.hospitalModel(),
                patientTumorLocation);

        CopyNumberAnalysis copyNumberAnalysis = analyzeCopyNumbers(purplePurityTsv, purpleGeneCnvTsv, patientTumorLocation);
        SomaticVariantAnalysis somaticVariantAnalysis =
                analyzeSomaticVariants(tumorSample, somaticVariantVcf, patientTumorLocation, copyNumberAnalysis.exomeGeneCopyNumbers());
        List<GermlineVariant> germlineVariantsToReport =
                analyzeGermlineVariants(tumorSample, bachelorCsv, copyNumberAnalysis, somaticVariantAnalysis);

        List<ReportableVariant> reportableVariants =
                ReportableVariantAnalyzer.mergeSomaticAndGermlineVariants(somaticVariantAnalysis.variantsToReport(),
                        somaticVariantAnalysis.driverCatalog(),
                        reportData.driverGeneView(),
                        germlineVariantsToReport,
                        reportData.germlineReportingModel(),
                        reportData.limsModel().germlineReportingChoice(tumorSample));

        SvAnalysis svAnalysis = analyzeStructuralVariants(linxFusionTsv, linxDisruptionTsv, copyNumberAnalysis, patientTumorLocation);
        ChordAnalysis chordAnalysis = analyzeChord(chordPredictionFile);

        String clinicalSummary = reportData.summaryModel().findSummaryForSample(tumorSample);

        List<EvidenceItem> allEvidenceItems = Lists.newArrayList();
        allEvidenceItems.addAll(somaticVariantAnalysis.evidenceItems());
        allEvidenceItems.addAll(copyNumberAnalysis.evidenceItems());
        allEvidenceItems.addAll(svAnalysis.evidenceItems());

        List<EvidenceItem> nonTrials = ReportableEvidenceItemFactory.extractNonTrials(allEvidenceItems);
        AnalysedPatientReport report = ImmutableAnalysedPatientReport.of(sampleReport,
                copyNumberAnalysis.hasReliablePurityFit(),
                copyNumberAnalysis.purity(),
                copyNumberAnalysis.ploidy(),
                clinicalSummary,
                nonTrials.stream().filter(EvidenceItem::isOnLabel).collect(Collectors.toList()),
                ClinicalTrialFactory.extractOnLabelTrials(allEvidenceItems),
                nonTrials.stream().filter(item -> !item.isOnLabel()).collect(Collectors.toList()),
                reportableVariants,
                somaticVariantAnalysis.microsatelliteIndelsPerMb(),
                somaticVariantAnalysis.tumorMutationalLoad(),
                somaticVariantAnalysis.tumorMutationalBurden(),
                chordAnalysis,
                copyNumberAnalysis.reportableGainsAndLosses(),
                svAnalysis.reportableFusions(),
                svAnalysis.reportableDisruptions(),
                circosFile,
                Optional.ofNullable(comments),
                Optional.ofNullable(correctTitle),
                reportData.signaturePath(),
                reportData.logoRVAPath(),
                reportData.logoCompanyPath());

        printReportState(report);

        return report;
    }

    @NotNull
    private CopyNumberAnalysis analyzeCopyNumbers(@NotNull String purplePurityTsv, @NotNull String purpleGeneCnvTsv,
            @Nullable PatientTumorLocation patientTumorLocation) throws IOException {
        PurityContext purityContext = FittedPurityFile.read(purplePurityTsv);

        DecimalFormat percentageFormat = new DecimalFormat("#'%'");
        LOGGER.info("Loaded purple sample data from {}", purplePurityTsv);
        LOGGER.info(" Purple purity: {}", percentageFormat.format(purityContext.bestFit().purity() * 100));
        LOGGER.info(" Purple average tumor ploidy: {}", purityContext.bestFit().ploidy());
        LOGGER.info(" Purple status: {}", purityContext.status());
        LOGGER.info(" WGD happened: {}", purityContext.wholeGenomeDuplication() ? "yes" : "no");

        List<GeneCopyNumber> exomeGeneCopyNumbers = GeneCopyNumberFile.read(purpleGeneCnvTsv);
        LOGGER.info("Loaded {} gene copy numbers from {}", exomeGeneCopyNumbers.size(), purpleGeneCnvTsv);

        return CopyNumberAnalyzer.run(purityContext, exomeGeneCopyNumbers, reportData.actionabilityAnalyzer(), patientTumorLocation);
    }

    @NotNull
    private SomaticVariantAnalysis analyzeSomaticVariants(@NotNull String sample, @NotNull String somaticVariantVcf,
            @Nullable PatientTumorLocation patientTumorLocation, @NotNull List<GeneCopyNumber> exomeGeneCopyNumbers) throws IOException {

        final List<SomaticVariant> variants = SomaticVariantFactory.passOnlyInstance().fromVCFFile(sample, somaticVariantVcf);
        LOGGER.info("Loaded {} PASS somatic variants from {}", variants.size(), somaticVariantVcf);

        return SomaticVariantAnalyzer.run(variants,
                reportData.driverGeneView(),
                reportData.actionabilityAnalyzer(),
                patientTumorLocation,
                exomeGeneCopyNumbers);
    }

    @NotNull
    private List<GermlineVariant> analyzeGermlineVariants(@NotNull String sample, @Nullable String bachelorCsv,
            @NotNull CopyNumberAnalysis copyNumberAnalysis, @NotNull SomaticVariantAnalysis somaticVariantAnalysis) throws IOException {
        if (bachelorCsv == null) {
            LOGGER.info("Skipping germline analysis - No bachelor CSV passed. Presumably no pathogenic germline variants found.");
            return Lists.newArrayList();
        }

        List<GermlineVariant> variants =
                BachelorFile.loadBachelorCsv(bachelorCsv).stream().filter(GermlineVariant::passFilter).collect(Collectors.toList());
        LOGGER.info("Loaded {} PASS germline variants from {}", variants.size(), bachelorCsv);

        LimsGermlineReportingChoice germlineChoice = reportData.limsModel().germlineReportingChoice(sample);
        if (germlineChoice == LimsGermlineReportingChoice.UNKNOWN) {
            LOGGER.info(" No germline reporting choice known. No germline variants will be reported!");
            return Lists.newArrayList();
        } else {
            LOGGER.info(" Patient has given the following germline consent: {}", germlineChoice);
            return FilterGermlineVariants.filterGermlineVariantsForReporting(variants,
                    reportData.driverGeneView(),
                    reportData.germlineReportingModel(),
                    copyNumberAnalysis.exomeGeneCopyNumbers(),
                    somaticVariantAnalysis.variantsToReport());
        }
    }

    @NotNull
    private SvAnalysis analyzeStructuralVariants(@NotNull String linxFusionTsv, @NotNull String linxDisruptionTsv,
            @NotNull CopyNumberAnalysis copyNumberAnalysis, @Nullable PatientTumorLocation patientTumorLocation) throws IOException {
        List<ReportableGeneFusion> fusions = ReportableGeneFusionFile.read(linxFusionTsv);
        LOGGER.info("Loaded {} fusions from {}", fusions.size(), linxFusionTsv);

        List<ReportableDisruption> disruptions = ReportableDisruptionFile.read(linxDisruptionTsv);
        LOGGER.info("Loaded {} disruptions from {}", disruptions.size(), linxDisruptionTsv);

        return SvAnalyzer.run(fusions,
                disruptions,
                copyNumberAnalysis.exomeGeneCopyNumbers(),
                reportData.actionabilityAnalyzer(),
                patientTumorLocation);
    }

    @NotNull
    private static ChordAnalysis analyzeChord(@NotNull String chordPredictionFile) throws IOException {
        ChordAnalysis chord = ChordFileReader.read(chordPredictionFile);
        LOGGER.info("Loaded CHORD analysis from {}", chordPredictionFile);
        return chord;
    }

    private static void printReportState(@NotNull AnalysedPatientReport report) {
        LocalDate tumorArrivalDate = report.sampleReport().tumorArrivalDate();
        String formattedTumorArrivalDate =
                tumorArrivalDate != null ? DateTimeFormatter.ofPattern("dd-MMM-yyyy").format(tumorArrivalDate) : "N/A";

        LOGGER.info("Printing clinical and laboratory data for {}", report.sampleReport().sampleId());
        LOGGER.info(" Tumor sample arrived at HMF on {}", formattedTumorArrivalDate);
        LOGGER.info(" Primary tumor location: {} ({})",
                report.sampleReport().primaryTumorLocationString(),
                report.sampleReport().cancerSubTypeString());
        LOGGER.info(" Shallow seq purity: {}", report.sampleReport().purityShallowSeq());
        LOGGER.info(" Lab SOPs used: {}", report.sampleReport().labProcedures());
        LOGGER.info(" Clinical summary present: {}", (!report.clinicalSummary().isEmpty() ? "yes" : "no"));

        List<ReportableVariant> variantsWithNotify =
                report.reportableVariants().stream().filter(ReportableVariant::notifyClinicalGeneticist).collect(Collectors.toList());
        LOGGER.info("Printing genomic analysis results for {}:", report.sampleReport().sampleId());
        LOGGER.info(" Somatic variants to report: {}", report.reportableVariants().size());
        LOGGER.info("  Variants for which to notify clinical geneticist: {}", variantsWithNotify.size());
        LOGGER.info(" Microsatellite Indels per Mb: {}", report.microsatelliteIndelsPerMb());
        LOGGER.info(" Tumor mutational load: {}", report.tumorMutationalLoad());
        LOGGER.info(" Tumor mutational burden: {}", report.tumorMutationalBurden());
        LOGGER.info(" CHORD analysis HRD prediction: {}", report.chordAnalysis().hrdValue());
        LOGGER.info(" Number of gains and losses to report: {}", report.gainsAndLosses().size());
        LOGGER.info(" Gene fusions to report : {}", report.geneFusions().size());
        LOGGER.info(" Gene disruptions to report : {}", report.geneDisruptions().size());

        LOGGER.info("Printing actionability results for {}", report.sampleReport().sampleId());
        LOGGER.info(" Tumor-specific evidence items found: {}", report.tumorSpecificEvidence().size());
        LOGGER.info(" Off-label evidence items found: {}", report.offLabelEvidence().size());
        LOGGER.info(" Clinical trials matched to molecular profile: {}", report.clinicalTrials().size());
    }
}
