package com.hartwig.hmftools.patientreporter;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocationFunctions;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingChoice;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.ClonalityCutoffKernel;
import com.hartwig.hmftools.common.variant.ClonalityFactory;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariantFactory;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariant;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariantFactory;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.patientreporter.actionability.ClinicalTrialFactory;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberAnalysis;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberAnalyzer;
import com.hartwig.hmftools.patientreporter.variants.germline.FilterGermlineVariants;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineVariant;
import com.hartwig.hmftools.patientreporter.structural.SvAnalysis;
import com.hartwig.hmftools.patientreporter.structural.SvAnalyzer;
import com.hartwig.hmftools.patientreporter.variants.ReportableVariant;
import com.hartwig.hmftools.patientreporter.variants.ReportableVariantAnalyzer;
import com.hartwig.hmftools.patientreporter.variants.somatic.SomaticVariantAnalysis;
import com.hartwig.hmftools.patientreporter.variants.somatic.SomaticVariantAnalyzer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

class PatientReporter {
    private static final Logger LOGGER = LogManager.getLogger(PatientReporter.class);

    @NotNull
    private final BaseReportData baseReportData;
    @NotNull
    private final SequencedReportData sequencedReportData;
    @NotNull
    private final SvAnalyzer svAnalyzerModel;

    PatientReporter(@NotNull final BaseReportData baseReportData, @NotNull final SequencedReportData sequencedReportData,
            @NotNull final SvAnalyzer svAnalyzerModel) {
        this.baseReportData = baseReportData;
        this.sequencedReportData = sequencedReportData;
        this.svAnalyzerModel = svAnalyzerModel;
    }

    @NotNull
    public AnalysedPatientReport run(@NotNull String runDirectory, @Nullable String comments) throws IOException {
        final RunContext run = ProductionRunContextFactory.fromRunDirectory(runDirectory);
        assert run.isSomaticRun();

        final String tumorSample = run.tumorSample();

        final PatientTumorLocation patientTumorLocation =
                PatientTumorLocationFunctions.findPatientTumorLocationForSample(baseReportData.patientTumorLocations(), tumorSample);

        final CopyNumberAnalysis copyNumberAnalysis = analyzeCopyNumbers(run, patientTumorLocation);
        final SomaticVariantAnalysis somaticVariantAnalysis = analyzeSomaticVariants(run, copyNumberAnalysis, patientTumorLocation);
        final List<GermlineVariant> germlineVariantsToReport = analyzeGermlineVariants(run, copyNumberAnalysis, somaticVariantAnalysis);

        final List<ReportableVariant> reportableVariants =
                ReportableVariantAnalyzer.mergeSomaticAndGermlineVariants(somaticVariantAnalysis.variantsToReport(),
                        somaticVariantAnalysis.driverCatalog(),
                        sequencedReportData.panelGeneModel().geneDriverCategoryMap(),
                        sequencedReportData.panelGeneModel().drupActionableGenes(),
                        germlineVariantsToReport,
                        sequencedReportData.germlineReportingModel(),
                        baseReportData.limsModel().germlineReportingChoice(tumorSample));

        final SvAnalysis svAnalysis = analyzeStructuralVariants(copyNumberAnalysis, patientTumorLocation, svAnalyzerModel);
        final ChordAnalysis chordAnalysis = analyzeChord(run);

        LOGGER.info("Printing analysis results:");
        LOGGER.info(" Somatic variants to report : " + somaticVariantAnalysis.variantsToReport().size());
        LOGGER.info(" Germline variants to report: " + germlineVariantsToReport.size());
        LOGGER.info("  Total number of reportable variants: " + reportableVariants.size());
        LOGGER.info(" Microsatellite Indels per Mb: " + somaticVariantAnalysis.microsatelliteIndelsPerMb());
        LOGGER.info(" Tumor mutational load: " + somaticVariantAnalysis.tumorMutationalLoad());
        LOGGER.info(" Tumor mutational burden: " + somaticVariantAnalysis.tumorMutationalBurden());
        LOGGER.info(" CHORD analysis HRD prediction: " + chordAnalysis.hrdValue());
        LOGGER.info(" Copy number events to report: " + copyNumberAnalysis.reportableGeneCopyNumbers().size());
        LOGGER.info(" Gene fusions to report : " + svAnalysis.reportableFusions().size());
        LOGGER.info(" Gene disruptions to report : " + svAnalysis.reportableDisruptions().size());
        LOGGER.info("Printing actionability results (including off-label trials):");
        LOGGER.info(" Evidence items found based on variants: " + somaticVariantAnalysis.evidenceItems().size());
        LOGGER.info(" Evidence items found based on copy numbers: " + copyNumberAnalysis.evidenceItems().size());
        LOGGER.info(" Evidence items found based on fusions: " + svAnalysis.evidenceItems().size());

        final List<EvidenceItem> allEvidenceItems = Lists.newArrayList();
        allEvidenceItems.addAll(somaticVariantAnalysis.evidenceItems());
        allEvidenceItems.addAll(copyNumberAnalysis.evidenceItems());
        allEvidenceItems.addAll(svAnalysis.evidenceItems());

        final SampleReport sampleReport = SampleReportFactory.fromLimsAndHospitalModel(tumorSample,
                run.refSample(),
                baseReportData.limsModel(),
                baseReportData.hospitalModel(),
                patientTumorLocation);

        final List<EvidenceItem> nonTrials = ReportableEvidenceItemFactory.extractNonTrials(allEvidenceItems);

        return ImmutableAnalysedPatientReport.of(sampleReport,
                copyNumberAnalysis.fittedPurity().purity(),
                copyNumberAnalysis.status() != FittedPurityStatus.NO_TUMOR,
                copyNumberAnalysis.fittedPurity().ploidy(),
                nonTrials.stream().filter(EvidenceItem::isOnLabel).collect(Collectors.toList()),
                ClinicalTrialFactory.extractOnLabelTrials(allEvidenceItems),
                nonTrials.stream().filter(item -> !item.isOnLabel()).collect(Collectors.toList()),
                reportableVariants,
                somaticVariantAnalysis.microsatelliteIndelsPerMb(),
                somaticVariantAnalysis.tumorMutationalLoad(),
                somaticVariantAnalysis.tumorMutationalBurden(),
                chordAnalysis,
                copyNumberAnalysis.reportableGeneCopyNumbers(),
                svAnalysis.reportableFusions(),
                svAnalysis.reportableDisruptions(),
                PatientReporterFileLoader.findCircosPlotPath(runDirectory, tumorSample),
                Optional.ofNullable(comments),
                baseReportData.signaturePath(),
                baseReportData.logoRVAPath());
    }

    @NotNull
    private CopyNumberAnalysis analyzeCopyNumbers(@NotNull RunContext run, @Nullable PatientTumorLocation patientTumorLocation)
            throws IOException {
        final String runDirectory = run.runDirectory();
        final String sample = run.tumorSample();

        LOGGER.info("Loading purple data for sample " + sample);
        final PurityContext purityContext = PatientReporterFileLoader.loadPurity(runDirectory, sample);
        LOGGER.info(" Purple purity " + purityContext.bestFit().purity());
        LOGGER.info(" Purple average tumor ploidy: " + purityContext.bestFit().ploidy());
        LOGGER.info(" Purple status " + purityContext.status());

        final List<PurpleCopyNumber> purpleCopyNumbers = PatientReporterFileLoader.loadPurpleCopyNumbers(runDirectory, sample);
        LOGGER.info(" " + purpleCopyNumbers.size() + " purple copy number regions loaded for sample " + sample);

        final List<GeneCopyNumber> exomeGeneCopyNumbers = PatientReporterFileLoader.loadPurpleGeneCopyNumbers(runDirectory, sample);

        LOGGER.info("Analyzing purple copy numbers");
        return CopyNumberAnalyzer.run(purityContext,
                purpleCopyNumbers,
                exomeGeneCopyNumbers,
                sequencedReportData.panelGeneModel(),
                sequencedReportData.actionabilityAnalyzer(),
                patientTumorLocation);
    }

    @NotNull
    private SomaticVariantAnalysis analyzeSomaticVariants(@NotNull RunContext run, @NotNull CopyNumberAnalysis copyNumberAnalysis,
            @Nullable PatientTumorLocation patientTumorLocation) throws IOException {
        final String runDirectory = run.runDirectory();
        final String sample = run.tumorSample();

        LOGGER.info("Loading somatic variants for sample " + sample);
        final List<SomaticVariant> variants =
                PatientReporterFileLoader.loadPassedSomaticVariants(runDirectory, sample, sequencedReportData.somaticVariantEnrichment());
        LOGGER.info(" " + variants.size() + " PASS somatic variants loaded for sample " + sample);

        LOGGER.info("Enriching somatic variants");
        final List<EnrichedSomaticVariant> enrichedSomaticVariants = enrich(variants,
                copyNumberAnalysis,
                sequencedReportData.highConfidenceRegions(),
                sequencedReportData.refGenomeFastaFile());

        LOGGER.info("Analyzing somatic variants");
        return SomaticVariantAnalyzer.run(enrichedSomaticVariants,
                sequencedReportData.panelGeneModel().somaticVariantGenes(),
                sequencedReportData.panelGeneModel().geneDriverCategoryMap(),
                sequencedReportData.actionabilityAnalyzer(),
                patientTumorLocation);
    }

    @NotNull
    private static List<EnrichedSomaticVariant> enrich(@NotNull List<SomaticVariant> variants,
            @NotNull CopyNumberAnalysis copyNumberAnalysis, @NotNull Multimap<String, GenomeRegion> highConfidenceRegions,
            @NotNull IndexedFastaSequenceFile refGenomeFastaFile) {
        final PurityAdjuster purityAdjuster = new PurityAdjuster(copyNumberAnalysis.gender(), copyNumberAnalysis.fittedPurity());
        final PurityAdjustedSomaticVariantFactory purityAdjustedFactory =
                new PurityAdjustedSomaticVariantFactory(purityAdjuster, copyNumberAnalysis.copyNumbers(), Collections.emptyList());
        final List<PurityAdjustedSomaticVariant> purityAdjustedSomaticVariants = purityAdjustedFactory.create(variants);

        final double clonalPloidy = ClonalityCutoffKernel.clonalCutoff(purityAdjustedSomaticVariants);
        final ClonalityFactory clonalityFactory = new ClonalityFactory(purityAdjuster, clonalPloidy);

        final EnrichedSomaticVariantFactory enrichedSomaticFactory =
                new EnrichedSomaticVariantFactory(highConfidenceRegions, refGenomeFastaFile, clonalityFactory);

        return enrichedSomaticFactory.enrich(purityAdjustedSomaticVariants);
    }

    @NotNull
    private List<GermlineVariant> analyzeGermlineVariants(@NotNull RunContext run, @NotNull CopyNumberAnalysis copyNumberAnalysis,
            @NotNull SomaticVariantAnalysis somaticVariantAnalysis) throws IOException {
        final String runDirectory = run.runDirectory();
        final String sample = run.tumorSample();

        LOGGER.info("Loading germline variants...");
        List<GermlineVariant> variants = PatientReporterFileLoader.loadPassedGermlineVariants(runDirectory, sample);
        LOGGER.info(" " + variants.size() + " PASS germline variants loaded for sample " + sample);

        LimsGermlineReportingChoice germlineChoice = baseReportData.limsModel().germlineReportingChoice(sample);
        if (germlineChoice == LimsGermlineReportingChoice.UNKNOWN) {
            LOGGER.info(" No germline reporting choice known. No germline variants will be reported!");
            return Lists.newArrayList();
        } else {
            return FilterGermlineVariants.filterGermlineVariantsForReporting(variants,
                    sequencedReportData.germlineReportingModel(),
                    sequencedReportData.panelGeneModel().geneDriverCategoryMap(),
                    copyNumberAnalysis.exomeGeneCopyNumbers(),
                    somaticVariantAnalysis.variantsToReport());
        }
    }

    @NotNull
    private SvAnalysis analyzeStructuralVariants(@NotNull CopyNumberAnalysis copyNumberAnalysis,
            @Nullable PatientTumorLocation patientTumorLocation, @NotNull SvAnalyzer svAnalyzer) {
        return svAnalyzer.run(sequencedReportData.panelGeneModel(),
                copyNumberAnalysis.exomeGeneCopyNumbers(),
                sequencedReportData.actionabilityAnalyzer(),
                patientTumorLocation);
    }

    @NotNull
    private static ChordAnalysis analyzeChord(@NotNull RunContext run) throws IOException {
        return PatientReporterFileLoader.loadChordFile(run.runDirectory(), run.tumorSample());
    }
}
