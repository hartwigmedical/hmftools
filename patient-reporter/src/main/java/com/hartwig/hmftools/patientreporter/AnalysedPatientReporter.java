package com.hartwig.hmftools.patientreporter;

import java.io.IOException;
import java.text.DecimalFormat;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
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
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableDisruption;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableDisruptionFile;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusionFile;
import com.hartwig.hmftools.common.variant.structural.linx.LinxViralInsertFile;
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
import com.hartwig.hmftools.patientreporter.viralInsertion.ImmutableViralInsertion;
import com.hartwig.hmftools.patientreporter.viralInsertion.ViralInsertion;

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
    AnalysedPatientReport run(@NotNull SampleMetadata sampleMetadata, @NotNull String purplePurityTsv, @NotNull String purpleGeneCnvTsv,
            @NotNull String somaticVariantVcf, @NotNull String linxFusionTsv, @NotNull String linxDisruptionTsv,
            @NotNull String bachelorCsv, @NotNull String chordPredictionFile, @NotNull String circosFile,
            @NotNull String linxViralInsertionFile, @Nullable String comments, boolean correctedReport) throws IOException {
        PatientTumorLocation patientTumorLocation =
                PatientTumorLocationFunctions.findPatientTumorLocationForSample(reportData.patientTumorLocations(),
                        sampleMetadata.tumorSampleId());

        SampleReport sampleReport = SampleReportFactory.fromLimsAndHospitalModel(sampleMetadata,
                reportData.limsModel(),
                reportData.hospitalModel(),
                patientTumorLocation);

        CopyNumberAnalysis copyNumberAnalysis = analyzeCopyNumbers(purplePurityTsv, purpleGeneCnvTsv, patientTumorLocation);
        SomaticVariantAnalysis somaticVariantAnalysis = analyzeSomaticVariants(sampleMetadata.tumorSampleId(),
                somaticVariantVcf,
                patientTumorLocation,
                copyNumberAnalysis.exomeGeneCopyNumbers());
        List<GermlineVariant> germlineVariantsToReport =
                analyzeGermlineVariants(sampleMetadata.tumorSampleBarcode(), bachelorCsv, copyNumberAnalysis, somaticVariantAnalysis);

        List<ReportableVariant> reportableVariants =
                ReportableVariantAnalyzer.mergeSomaticAndGermlineVariants(somaticVariantAnalysis.variantsToReport(),
                        somaticVariantAnalysis.driverCatalog(),
                        reportData.driverGeneView(),
                        germlineVariantsToReport,
                        reportData.germlineReportingModel(),
                        reportData.limsModel().germlineReportingChoice(sampleMetadata.tumorSampleBarcode()));

        SvAnalysis svAnalysis = analyzeStructuralVariants(linxFusionTsv, linxDisruptionTsv, copyNumberAnalysis, patientTumorLocation);
        List<ViralInsertion> viralInsertions = analyzeViralInsertions(linxViralInsertionFile);
        ChordAnalysis chordAnalysis = analyzeChord(chordPredictionFile);

        String clinicalSummary = reportData.summaryModel().findSummaryForSample(sampleMetadata.tumorSampleId());

        List<EvidenceItem> allEvidenceItems = Lists.newArrayList();
        allEvidenceItems.addAll(somaticVariantAnalysis.evidenceItems());
        allEvidenceItems.addAll(copyNumberAnalysis.evidenceItems());
        allEvidenceItems.addAll(svAnalysis.evidenceItems());

        List<EvidenceItem> nonTrials = ReportableEvidenceItemFactory.extractNonTrials(allEvidenceItems);
        AnalysedPatientReport report = ImmutableAnalysedPatientReport.builder()
                .sampleReport(sampleReport)
                .hasReliablePurityFit(copyNumberAnalysis.hasReliablePurityFit())
                .impliedPurity(copyNumberAnalysis.purity())
                .averageTumorPloidy(copyNumberAnalysis.ploidy())
                .clinicalSummary(clinicalSummary)
                .tumorSpecificEvidence(nonTrials.stream().filter(EvidenceItem::isOnLabel).collect(Collectors.toList()))
                .clinicalTrials(ClinicalTrialFactory.extractOnLabelTrials(allEvidenceItems))
                .offLabelEvidence(nonTrials.stream().filter(item -> !item.isOnLabel()).collect(Collectors.toList()))
                .reportableVariants(reportableVariants)
                .microsatelliteIndelsPerMb(somaticVariantAnalysis.microsatelliteIndelsPerMb())
                .tumorMutationalLoad(somaticVariantAnalysis.tumorMutationalLoad())
                .tumorMutationalBurden(somaticVariantAnalysis.tumorMutationalBurden())
                .chordAnalysis(chordAnalysis)
                .gainsAndLosses(copyNumberAnalysis.reportableGainsAndLosses())
                .geneFusions(svAnalysis.reportableFusions())
                .geneDisruptions(svAnalysis.reportableDisruptions())
                .viralInsertion(viralInsertions)
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
    private List<ViralInsertion> analyzeViralInsertions(@NotNull String linxViralInsertionFile) throws IOException {
        List<LinxViralInsertFile> viralInsertFileList = LinxViralInsertFile.read(linxViralInsertionFile);
        LOGGER.info("Loaded {} viral insertions from {}", viralInsertFileList.size(), linxViralInsertionFile);

        Map<AnalysedPatientReporter.VirusKey, List<LinxViralInsertFile>> itemsPerKey = Maps.newHashMap();
        for (LinxViralInsertFile viralInsertion : viralInsertFileList) {
            AnalysedPatientReporter.VirusKey key = new AnalysedPatientReporter.VirusKey(viralInsertion.VirusId);
            List<LinxViralInsertFile> items = itemsPerKey.get(key);

            if (items == null) {
                items = Lists.newArrayList();
            }
            items.add(viralInsertion);
            itemsPerKey.put(key, items);
        }

        List<ViralInsertion> viralInsertions = Lists.newArrayList();
        String virusName = Strings.EMPTY;
        int count = 0;
        for (Map.Entry<AnalysedPatientReporter.VirusKey, List<LinxViralInsertFile>> entry : itemsPerKey.entrySet()) {
            List<LinxViralInsertFile> itemsForKey = entry.getValue();
            for (LinxViralInsertFile virus: itemsForKey) {
                 virusName = virus.VirusName;
            }

             count = itemsForKey.size();

            viralInsertions.add(ImmutableViralInsertion.builder().virus(virusName).countVirus(Integer.toString(count)).build());
        }
        return viralInsertions;
    }

    private static class VirusKey {

        @NotNull
        private final String virusId;

        private VirusKey(@NotNull final String virusId) {
            this.virusId = virusId;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) {
                return true;
            }
            if (o == null || getClass() != o.getClass()) {
                return false;
            }
            final AnalysedPatientReporter.VirusKey key = (AnalysedPatientReporter.VirusKey) o;
            return Objects.equals(virusId, key.virusId);
        }

        @Override
        public int hashCode() {
            return Objects.hash(virusId);
        }
    }

    @NotNull
    private CopyNumberAnalysis analyzeCopyNumbers(@NotNull String purplePurityTsv, @NotNull String purpleGeneCnvTsv,
            @Nullable PatientTumorLocation patientTumorLocation) throws IOException {
        PurityContext purityContext = FittedPurityFile.read(purplePurityTsv);

        LOGGER.info("Loaded purple sample data from {}", purplePurityTsv);
        LOGGER.info(" Purple purity: {}", new DecimalFormat("#'%'").format(purityContext.bestFit().purity() * 100));
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
    private List<GermlineVariant> analyzeGermlineVariants(@NotNull String sampleBarcode, @NotNull String bachelorCsv,
            @NotNull CopyNumberAnalysis copyNumberAnalysis, @NotNull SomaticVariantAnalysis somaticVariantAnalysis) throws IOException {

        List<GermlineVariant> variants =
                BachelorFile.loadBachelorCsv(bachelorCsv).stream().filter(GermlineVariant::passFilter).collect(Collectors.toList());
        LOGGER.info("Loaded {} PASS germline variants from {}", variants.size(), bachelorCsv);

        LimsGermlineReportingChoice germlineChoice = reportData.limsModel().germlineReportingChoice(sampleBarcode);
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
