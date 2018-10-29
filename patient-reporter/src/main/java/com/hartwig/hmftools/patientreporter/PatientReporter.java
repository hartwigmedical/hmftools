package com.hartwig.hmftools.patientreporter;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeMappingReading;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.collect.Multimaps;
import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
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
import com.hartwig.hmftools.common.variant.enrich.SomaticEnrichment;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariantFactory;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFileLoader;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.patientreporter.actionability.ActionabilityVariantAnalyzer;
import com.hartwig.hmftools.patientreporter.actionability.ClinicalTrialFactory;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;
import com.hartwig.hmftools.patientreporter.chord.ChordAnalysis;
import com.hartwig.hmftools.patientreporter.copynumber.ImmutablePurpleAnalysis;
import com.hartwig.hmftools.patientreporter.copynumber.PurpleAnalysis;
import com.hartwig.hmftools.patientreporter.genepanel.GeneModel;
import com.hartwig.hmftools.patientreporter.germline.GermlineVariant;
import com.hartwig.hmftools.patientreporter.structural.ImmutableReportableStructuralVariantAnalysis;
import com.hartwig.hmftools.patientreporter.structural.ReportableGeneDisruption;
import com.hartwig.hmftools.patientreporter.structural.ReportableGeneDisruptionFactory;
import com.hartwig.hmftools.patientreporter.structural.ReportableGeneFusion;
import com.hartwig.hmftools.patientreporter.structural.ReportableGeneFusionFactory;
import com.hartwig.hmftools.patientreporter.structural.ReportableStructuralVariantAnalysis;
import com.hartwig.hmftools.patientreporter.variants.SomaticVariantAnalysis;
import com.hartwig.hmftools.patientreporter.variants.SomaticVariantAnalyzer;
import com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalysis;
import com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalyzer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
abstract class PatientReporter {
    private static final Logger LOGGER = LogManager.getLogger(PatientReporter.class);

    @NotNull
    public abstract BaseReportData baseReportData();

    @NotNull
    public abstract SequencedReportData sequencedReportData();

    @NotNull
    public abstract StructuralVariantAnalyzer structuralVariantAnalyzer();

    @NotNull
    public AnalysedPatientReport run(@NotNull String runDirectory, boolean doReportGermline, @Nullable String comments) throws IOException {
        final RunContext run = ProductionRunContextFactory.fromRunDirectory(runDirectory);
        assert run.isSomaticRun();

        final String tumorSample = run.tumorSample();
        final PatientTumorLocation patientTumorLocation =
                PatientReporterFileLoader.extractPatientTumorLocation(baseReportData().patientTumorLocations(), tumorSample);

        final PurpleAnalysis purpleAnalysis = analyzePurpleCopyNumbers(run,
                sequencedReportData().actionabilityAnalyzer(),
                patientTumorLocation,
                sequencedReportData().panelGeneModel());
        final List<GeneCopyNumber> reportableGeneCopynumbers =
                purpleAnalysis.reportableGeneCopyNumbers(sequencedReportData().panelGeneModel());

        final ReportableStructuralVariantAnalysis reportableStructuralVariantAnalysis = analyzeStructuralVariants(run,
                purpleAnalysis,
                structuralVariantAnalyzer(),
                patientTumorLocation,
                sequencedReportData().actionabilityAnalyzer());

        final ChordAnalysis chordAnalysis = analyzeChord(run);

        final SomaticVariantAnalysis somaticVariantAnalysis = analyzeSomaticVariants(run,
                purpleAnalysis,
                sequencedReportData().somaticVariantEnrichment(),
                sequencedReportData().panelGeneModel(),
                sequencedReportData().highConfidenceRegions(),
                sequencedReportData().refGenomeFastaFile(),
                patientTumorLocation,
                sequencedReportData().actionabilityAnalyzer());

        final List<GermlineVariant> germlineVariants = doReportGermline ? analyzeGermlineVariants(run) : null;

        LOGGER.info("Printing analysis results:");
        LOGGER.info(
                " Number of somatic variants to report : " + Integer.toString(somaticVariantAnalysis.reportableSomaticVariants().size()));
        LOGGER.info(" Microsatellite analysis results: " + Double.toString(somaticVariantAnalysis.microsatelliteIndelsPerMb())
                + " indels per MB");
        LOGGER.info(" Tumor mutational load: " + Integer.toString(somaticVariantAnalysis.tumorMutationalLoad()));
        LOGGER.info(" Tumor mutational burden: " + Double.toString(somaticVariantAnalysis.tumorMutationalBurden()) + " mutations per MB");
        LOGGER.info(" CHORD analysis HRD prediction: " + Double.toString(chordAnalysis.hrdValue()));
        LOGGER.info(" Number of germline variants to report : " + Integer.toString(germlineVariants != null ? germlineVariants.size() : 0));
        LOGGER.info(" Number of copy number events to report: " + Integer.toString(reportableGeneCopynumbers.size()));
        LOGGER.info(
                " Number of gene fusions to report : " + Integer.toString(reportableStructuralVariantAnalysis.reportableFusions().size()));
        LOGGER.info(
                " Number of gene disruptions to report : " + Integer.toString(reportableStructuralVariantAnalysis.reportableDisruptions()
                        .size()));

        final SampleReport sampleReport = ImmutableSampleReport.of(tumorSample,
                patientTumorLocation,
                baseReportData().limsModel().tumorPercentageForSample(tumorSample),
                baseReportData().limsModel().arrivalDateForSample(tumorSample),
                baseReportData().limsModel().arrivalDateForSample(run.refSample()),
                baseReportData().limsModel().labProceduresForSample(tumorSample),
                baseReportData().centerModel().getAddresseeStringForSample(tumorSample));

        List<EvidenceItem> evidenceItems = Lists.newArrayList();
        for (Map.Entry<EnrichedSomaticVariant, List<EvidenceItem>> evidencePerVariant : somaticVariantAnalysis.evidencePerVariant()
                .entrySet()) {
            evidenceItems.addAll(evidencePerVariant.getValue());
        }

        for (Map.Entry<GeneCopyNumber, List<EvidenceItem>> evidencePerCNV : purpleAnalysis.evidencePerGeneCopyNumber().entrySet()) {
            evidenceItems.addAll(evidencePerCNV.getValue());
        }

        for (Map.Entry<GeneFusion, List<EvidenceItem>> evidencePerFusion : reportableStructuralVariantAnalysis.evidencePerFusion()
                .entrySet()) {
            evidenceItems.addAll(evidencePerFusion.getValue());
        }

        LOGGER.info("Printing actionability results:");
        LOGGER.info(" Number of evidence items based on variants: " + Integer.toString(somaticVariantAnalysis.evidencePerVariant().size()));
        LOGGER.info(
                " Number of evidence items based on copy numbers: " + Integer.toString(purpleAnalysis.evidencePerGeneCopyNumber().size()));
        LOGGER.info(
                " Number of evidence items based on fusions: " + Integer.toString(reportableStructuralVariantAnalysis.evidencePerFusion()
                        .size()));

        return ImmutableAnalysedPatientReport.of(sampleReport,
                purpleAnalysis.fittedPurity().purity(),
                purpleAnalysis.status() != FittedPurityStatus.NO_TUMOR,
                ReportableEvidenceItemFactory.filterEvidenceItemsForReporting(evidenceItems),
                ClinicalTrialFactory.extractTrials(evidenceItems),
                somaticVariantAnalysis.reportableSomaticVariants(),
                somaticVariantAnalysis.microsatelliteIndelsPerMb(),
                somaticVariantAnalysis.tumorMutationalLoad(),
                somaticVariantAnalysis.tumorMutationalBurden(),
                chordAnalysis,
                germlineVariants != null,
                germlineVariants != null ? germlineVariants : Lists.newArrayList(),
                reportableGeneCopynumbers,
                reportableStructuralVariantAnalysis.reportableFusions(),
                reportableStructuralVariantAnalysis.reportableDisruptions(),
                PatientReporterFileLoader.findCircosPlotPath(runDirectory, tumorSample),
                Optional.ofNullable(comments),
                baseReportData().signaturePath(),
                baseReportData().logoPath());
    }

    @NotNull
    private static PurpleAnalysis analyzePurpleCopyNumbers(@NotNull RunContext run, @NotNull ActionabilityAnalyzer actionabilityAnalyzer,
            @Nullable PatientTumorLocation patientTumorLocation, @NotNull GeneModel panelGeneModel) throws IOException {
        final String runDirectory = run.runDirectory();
        final String sample = run.tumorSample();

        LOGGER.info("Loading purple data for sample " + sample);
        final PurityContext purityContext = PatientReporterFileLoader.loadPurity(runDirectory, sample);

        final List<PurpleCopyNumber> purpleCopyNumbers = PatientReporterFileLoader.loadPurpleCopyNumbers(runDirectory, sample);
        LOGGER.info(" " + purpleCopyNumbers.size() + " purple copy number regions loaded for sample " + sample);

        final List<GeneCopyNumber> geneCopyNumbers = PatientReporterFileLoader.loadPurpleGeneCopyNumbers(runDirectory, sample);

        LOGGER.info("Determining evidence for copy numbers");
        Set<String> actionableGenesCNVS = actionabilityAnalyzer.cnvAnalyzer().actionableGenes();
        final String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : Strings.EMPTY;
        CancerTypeMappingReading cancerTypeMappingReading = CancerTypeMappingReading.readingFile();
        String doidsPrimaryTumorLocation = cancerTypeMappingReading.doidsForPrimaryTumorLocation(primaryTumorLocation);

        final Map<GeneCopyNumber, List<EvidenceItem>> evidencePerGeneCopyNumber = ActionabilityVariantAnalyzer.findEvidenceForCopyNumber(
                actionableGenesCNVS,
                geneCopyNumbers,
                doidsPrimaryTumorLocation,
                actionabilityAnalyzer,
                purityContext.bestFit().ploidy());

        final List<GeneCopyNumber> panelGeneCopyNumbers = geneCopyNumbers.stream()
                .filter(geneCopyNumber -> panelGeneModel.cnvGenePanel().contains(geneCopyNumber.gene()))
                .collect(Collectors.toList());

        return ImmutablePurpleAnalysis.builder()
                .gender(purityContext.gender())
                .status(purityContext.status())
                .fittedPurity(purityContext.bestFit())
                .fittedScorePurity(purityContext.score())
                .copyNumbers(purpleCopyNumbers)
                .geneCopyNumbers(geneCopyNumbers)
                .evidencePerGeneCopyNumber(evidencePerGeneCopyNumber)
                .panelGeneCopyNumbers(panelGeneCopyNumbers)
                .build();
    }

    @NotNull
    private static SomaticVariantAnalysis analyzeSomaticVariants(@NotNull RunContext run, @NotNull PurpleAnalysis purpleAnalysis,
            @NotNull SomaticEnrichment somaticEnrichment, @NotNull GeneModel geneModel,
            @NotNull Multimap<String, GenomeRegion> highConfidenceRegions, @NotNull IndexedFastaSequenceFile refGenomeFastaFile,
            @Nullable PatientTumorLocation patientTumorLocation, @NotNull ActionabilityAnalyzer actionabilityAnalyzerData)
            throws IOException {
        final String runDirectory = run.runDirectory();
        final String sample = run.tumorSample();

        LOGGER.info("Loading somatic variants...");
        final List<SomaticVariant> variants = PatientReporterFileLoader.loadPassedSomaticVariants(runDirectory, sample, somaticEnrichment);
        LOGGER.info(" " + variants.size() + " PASS somatic variants loaded for sample " + sample);

        LOGGER.info("Enriching somatic variants");
        final List<EnrichedSomaticVariant> enrichedSomaticVariants =
                enrich(variants, purpleAnalysis, highConfidenceRegions, refGenomeFastaFile);

        LOGGER.info("Analyzing somatic variants....");
        return SomaticVariantAnalyzer.run(enrichedSomaticVariants,
                geneModel.somaticVariantGenePanel(),
                geneModel.geneDriverCategoryMap(),
                geneModel.drupActionableGenes().keySet(),
                patientTumorLocation,
                actionabilityAnalyzerData);
    }

    @NotNull
    private static List<EnrichedSomaticVariant> enrich(@NotNull List<SomaticVariant> variants, @NotNull PurpleAnalysis purpleAnalysis,
            @NotNull Multimap<String, GenomeRegion> highConfidenceRegions, @NotNull IndexedFastaSequenceFile refGenomeFastaFile) {
        final PurityAdjuster purityAdjuster = new PurityAdjuster(purpleAnalysis.gender(), purpleAnalysis.fittedPurity());
        final PurityAdjustedSomaticVariantFactory purityAdjustedFactory =
                new PurityAdjustedSomaticVariantFactory(purityAdjuster, purpleAnalysis.copyNumbers(), Collections.emptyList());
        final List<PurityAdjustedSomaticVariant> purityAdjustedSomaticVariants = purityAdjustedFactory.create(variants);

        final double clonalPloidy = ClonalityCutoffKernel.clonalCutoff(purityAdjustedSomaticVariants);
        final ClonalityFactory clonalityFactory = new ClonalityFactory(purityAdjuster, clonalPloidy);

        final EnrichedSomaticVariantFactory enrichedSomaticFactory =
                new EnrichedSomaticVariantFactory(highConfidenceRegions, refGenomeFastaFile, clonalityFactory);

        return enrichedSomaticFactory.enrich(purityAdjustedSomaticVariants);
    }

    @NotNull
    private static ReportableStructuralVariantAnalysis analyzeStructuralVariants(@NotNull RunContext run,
            @NotNull PurpleAnalysis purpleAnalysis, @NotNull StructuralVariantAnalyzer structuralVariantAnalyzer,
            @Nullable PatientTumorLocation patientTumorLocation, @NotNull ActionabilityAnalyzer actionabilityAnalyzer) throws IOException {
        final Path structuralVariantVCF = PatientReporterFileLoader.findStructuralVariantVCF(run.runDirectory());
        LOGGER.info("Loading structural variants...");
        final List<StructuralVariant> structuralVariants = StructuralVariantFileLoader.fromFile(structuralVariantVCF.toString(), true);

        LOGGER.info("Enriching structural variants with purple data.");
        final PurityAdjuster purityAdjuster = new PurityAdjuster(purpleAnalysis.gender(), purpleAnalysis.fittedPurity());
        final Multimap<Chromosome, PurpleCopyNumber> copyNumberMap = Multimaps.fromRegions(purpleAnalysis.copyNumbers());

        final List<EnrichedStructuralVariant> enrichedStructuralVariants =
                EnrichedStructuralVariantFactory.enrich(structuralVariants, purityAdjuster, copyNumberMap);

        LOGGER.info("Analysing structural variants...");
        final StructuralVariantAnalysis structuralVariantAnalysis = structuralVariantAnalyzer.run(enrichedStructuralVariants);

        final List<ReportableGeneFusion> reportableFusions =
                ReportableGeneFusionFactory.toReportableGeneFusions(structuralVariantAnalysis.reportableFusions());
        final List<ReportableGeneDisruption> reportableDisruptions = ReportableGeneDisruptionFactory.toReportableGeneDisruptions(
                structuralVariantAnalysis.reportableDisruptions(),
                purpleAnalysis.geneCopyNumbers());

        final String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : Strings.EMPTY;
        CancerTypeMappingReading cancerTypeMappingReading = CancerTypeMappingReading.readingFile();
        String doidsPrimaryTumorLocation = cancerTypeMappingReading.doidsForPrimaryTumorLocation(primaryTumorLocation);

        Set<String> actionableFusions = actionabilityAnalyzer.fusionAnalyzer().actionableGenes();

        LOGGER.info("Analyzing fusions for actionability");
        Map<GeneFusion, List<EvidenceItem>> evidencePerFusion = ActionabilityVariantAnalyzer.findEvidenceForFusions(actionableFusions,
                structuralVariantAnalysis.fusions(),
                doidsPrimaryTumorLocation,
                actionabilityAnalyzer);

        return ImmutableReportableStructuralVariantAnalysis.of(reportableFusions, reportableDisruptions, evidencePerFusion);
    }

    @Nullable
    private static List<GermlineVariant> analyzeGermlineVariants(@NotNull RunContext run) throws IOException {
        final String runDirectory = run.runDirectory();
        final String sample = run.tumorSample();

        LOGGER.info("Loading germline variants...");
        final List<GermlineVariant> variants = PatientReporterFileLoader.loadPassedGermlineVariants(runDirectory, sample);
        if (variants == null) {
            LOGGER.warn(" Could not load germline variants. Probably bachelor hasn't been run yet!");
        } else {
            LOGGER.info(" " + variants.size() + " PASS germline variants loaded for sample " + sample);
        }

        return variants;
    }

    @NotNull
    private static ChordAnalysis analyzeChord(@NotNull RunContext run) throws IOException {
        return PatientReporterFileLoader.loadChordFile(run.runDirectory(), run.tumorSample());
    }
}
