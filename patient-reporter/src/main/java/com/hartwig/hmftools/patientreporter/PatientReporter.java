package com.hartwig.hmftools.patientreporter;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Collections;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.collect.Multimaps;
import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocationFunctions;
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
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariantFactory;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFileLoader;
import com.hartwig.hmftools.patientreporter.actionability.ClinicalTrialFactory;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;
import com.hartwig.hmftools.patientreporter.chord.ChordAnalysis;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberAnalysis;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberAnalyzer;
import com.hartwig.hmftools.patientreporter.germline.GermlineVariant;
import com.hartwig.hmftools.patientreporter.structural.FusionDisruptionAnalysis;
import com.hartwig.hmftools.patientreporter.structural.FusionDisruptionAnalyzer;
import com.hartwig.hmftools.patientreporter.variants.SomaticVariantAnalysis;
import com.hartwig.hmftools.patientreporter.variants.SomaticVariantAnalyzer;
import com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalysis;
import com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalyzer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
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
                PatientTumorLocationFunctions.findPatientTumorLocationForSample(baseReportData().patientTumorLocations(), tumorSample);

        final CopyNumberAnalysis copyNumberAnalysis = analyzeCopyNumbers(run, patientTumorLocation);
        final SomaticVariantAnalysis somaticVariantAnalysis = analyzeSomaticVariants(run, copyNumberAnalysis, patientTumorLocation);
        final FusionDisruptionAnalysis fusionDisruptionAnalysis = analyzeStructuralVariants(run, copyNumberAnalysis, patientTumorLocation);
        final List<GermlineVariant> germlineVariants = doReportGermline ? analyzeGermlineVariants(run) : null;
        final ChordAnalysis chordAnalysis = analyzeChord(run);

        LOGGER.info("Printing analysis results:");
        LOGGER.info(" Somatic variants to report : " + somaticVariantAnalysis.reportableSomaticVariants().size());
        LOGGER.info(" Microsatellite Indels per Mb: " + somaticVariantAnalysis.microsatelliteIndelsPerMb());
        LOGGER.info(" Tumor mutational load: " + somaticVariantAnalysis.tumorMutationalLoad());
        LOGGER.info(" Tumor mutational burden: " + somaticVariantAnalysis.tumorMutationalBurden());
        LOGGER.info(" CHORD analysis HRD prediction: " + chordAnalysis.hrdValue());
        LOGGER.info(" Germline variants to report : " + Integer.toString(germlineVariants != null ? germlineVariants.size() : 0));
        LOGGER.info(" Copy number events to report: " + copyNumberAnalysis.reportableGeneCopyNumbers().size());
        LOGGER.info(" Gene fusions to report : " + fusionDisruptionAnalysis.reportableFusions().size());
        LOGGER.info(" Gene disruptions to report : " + fusionDisruptionAnalysis.reportableDisruptions().size());
        LOGGER.info("Printing actionability results:");
        LOGGER.info(" Evidence items found based on variants: " + somaticVariantAnalysis.evidenceItems().size());
        LOGGER.info(" Evidence items found based on copy numbers: " + copyNumberAnalysis.evidenceItems().size());
        LOGGER.info(" Evidence items found based on fusions: " + fusionDisruptionAnalysis.evidenceItems().size());

        final List<EvidenceItem> allEvidenceItems = Lists.newArrayList();
        allEvidenceItems.addAll(somaticVariantAnalysis.evidenceItems());
        allEvidenceItems.addAll(copyNumberAnalysis.evidenceItems());
        allEvidenceItems.addAll(fusionDisruptionAnalysis.evidenceItems());

        final SampleReport sampleReport = ImmutableSampleReport.of(tumorSample,
                patientTumorLocation,
                baseReportData().limsModel().tumorPercentageForSample(tumorSample),
                baseReportData().limsModel().arrivalDateForSample(tumorSample),
                baseReportData().limsModel().arrivalDateForSample(run.refSample()),
                baseReportData().limsModel().labProceduresForSample(tumorSample),
                baseReportData().centerModel().getAddresseeStringForSample(tumorSample));

        return ImmutableAnalysedPatientReport.of(sampleReport,
                copyNumberAnalysis.fittedPurity().purity(),
                copyNumberAnalysis.status() != FittedPurityStatus.NO_TUMOR,
                copyNumberAnalysis.fittedPurity().ploidy(),
                ReportableEvidenceItemFactory.extractClinicalEvidence(allEvidenceItems),
                ClinicalTrialFactory.extractTrials(allEvidenceItems),
                somaticVariantAnalysis.reportableSomaticVariants(),
                somaticVariantAnalysis.microsatelliteIndelsPerMb(),
                somaticVariantAnalysis.tumorMutationalLoad(),
                somaticVariantAnalysis.tumorMutationalBurden(),
                chordAnalysis,
                germlineVariants != null,
                germlineVariants != null ? germlineVariants : Lists.newArrayList(),
                copyNumberAnalysis.reportableGeneCopyNumbers(),
                fusionDisruptionAnalysis.reportableFusions(),
                fusionDisruptionAnalysis.reportableDisruptions(),
                PatientReporterFileLoader.findCircosPlotPath(runDirectory, tumorSample),
                Optional.ofNullable(comments),
                baseReportData().signaturePath(),
                baseReportData().logoRVAPath());
    }

    @NotNull
    private CopyNumberAnalysis analyzeCopyNumbers(@NotNull RunContext run, @Nullable PatientTumorLocation patientTumorLocation)
            throws IOException {
        final String runDirectory = run.runDirectory();
        final String sample = run.tumorSample();

        LOGGER.info("Loading purple data for sample " + sample);
        final PurityContext purityContext = PatientReporterFileLoader.loadPurity(runDirectory, sample);

        final List<PurpleCopyNumber> purpleCopyNumbers = PatientReporterFileLoader.loadPurpleCopyNumbers(runDirectory, sample);
        LOGGER.info(" " + purpleCopyNumbers.size() + " purple copy number regions loaded for sample " + sample);

        final List<GeneCopyNumber> exomeGeneCopyNumbers = PatientReporterFileLoader.loadPurpleGeneCopyNumbers(runDirectory, sample);

        LOGGER.info("Analyzing purple copy numbers");
        return CopyNumberAnalyzer.run(purityContext,
                purpleCopyNumbers,
                exomeGeneCopyNumbers,
                sequencedReportData().panelGeneModel(),
                sequencedReportData().actionabilityAnalyzer(),
                patientTumorLocation);
    }

    @NotNull
    private SomaticVariantAnalysis analyzeSomaticVariants(@NotNull RunContext run, @NotNull CopyNumberAnalysis copyNumberAnalysis,
            @Nullable PatientTumorLocation patientTumorLocation) throws IOException {
        final String runDirectory = run.runDirectory();
        final String sample = run.tumorSample();

        LOGGER.info("Loading somatic variants for sample " + sample);
        final List<SomaticVariant> variants =
                PatientReporterFileLoader.loadPassedSomaticVariants(runDirectory, sample, sequencedReportData().somaticVariantEnrichment());
        LOGGER.info(" " + variants.size() + " PASS somatic variants loaded for sample " + sample);

        LOGGER.info("Enriching somatic variants");
        final List<EnrichedSomaticVariant> enrichedSomaticVariants = enrich(variants,
                copyNumberAnalysis,
                sequencedReportData().highConfidenceRegions(),
                sequencedReportData().refGenomeFastaFile());

        LOGGER.info("Analyzing somatic variants");
        return SomaticVariantAnalyzer.run(enrichedSomaticVariants,
                sequencedReportData().panelGeneModel().somaticVariantGenePanel(),
                sequencedReportData().panelGeneModel().geneDriverCategoryMap(),
                sequencedReportData().panelGeneModel().drupActionableGenes().keySet(),
                sequencedReportData().actionabilityAnalyzer(),
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
    private FusionDisruptionAnalysis analyzeStructuralVariants(@NotNull RunContext run, @NotNull CopyNumberAnalysis copyNumberAnalysis,
            @Nullable PatientTumorLocation patientTumorLocation) throws IOException {
        final Path structuralVariantVCF = PatientReporterFileLoader.findStructuralVariantVCF(run.runDirectory());
        LOGGER.info("Loading structural variants...");
        final List<StructuralVariant> structuralVariants = StructuralVariantFileLoader.fromFile(structuralVariantVCF.toString(), true);
        LOGGER.info(" " + structuralVariants.size() + " structural variants loaded");

        LOGGER.info("Enriching structural variants with purple data");
        final PurityAdjuster purityAdjuster = new PurityAdjuster(copyNumberAnalysis.gender(), copyNumberAnalysis.fittedPurity());
        final Multimap<Chromosome, PurpleCopyNumber> copyNumberMap = Multimaps.fromRegions(copyNumberAnalysis.copyNumbers());

        final List<EnrichedStructuralVariant> enrichedStructuralVariants =
                new EnrichedStructuralVariantFactory(purityAdjuster, copyNumberMap).enrich(structuralVariants);

        LOGGER.info("Analyzing structural variants...");
        final StructuralVariantAnalysis structuralVariantAnalysis = structuralVariantAnalyzer().run(enrichedStructuralVariants);

        return FusionDisruptionAnalyzer.run(structuralVariantAnalysis,
                copyNumberAnalysis.exomeGeneCopyNumbers(),
                sequencedReportData().actionabilityAnalyzer(),
                patientTumorLocation);
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
