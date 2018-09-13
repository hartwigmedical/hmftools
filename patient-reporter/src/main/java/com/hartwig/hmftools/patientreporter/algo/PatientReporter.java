package com.hartwig.hmftools.patientreporter.algo;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.common.variant.ClonalityCutoffKernel;
import com.hartwig.hmftools.common.variant.ClonalityFactory;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariantFactory;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariant;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariantFactory;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFileLoader;
import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.BaseReporterData;
import com.hartwig.hmftools.patientreporter.HmfReporterData;
import com.hartwig.hmftools.patientreporter.ImmutableAnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.ImmutableSampleReport;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.copynumber.ImmutablePurpleAnalysis;
import com.hartwig.hmftools.patientreporter.copynumber.PurpleAnalysis;
import com.hartwig.hmftools.patientreporter.report.data.GeneDisruptionData;
import com.hartwig.hmftools.patientreporter.report.data.GeneFusionData;
import com.hartwig.hmftools.patientreporter.variants.SomaticVariantAnalysis;
import com.hartwig.hmftools.patientreporter.variants.SomaticVariantAnalyzer;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;
import com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalysis;
import com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalyzer;
import com.hartwig.hmftools.svannotation.annotations.GeneDisruption;
import com.hartwig.hmftools.svannotation.annotations.GeneFusion;
import com.hartwig.hmftools.svannotation.annotations.Transcript;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class PatientReporter {
    private static final Logger LOGGER = LogManager.getLogger(PatientReporter.class);

    @NotNull
    public abstract BaseReporterData baseReporterData();

    @NotNull
    public abstract HmfReporterData reporterData();

    @NotNull
    public abstract SomaticVariantAnalyzer variantAnalyzer();

    @NotNull
    public abstract StructuralVariantAnalyzer structuralVariantAnalyzer();

    @NotNull
    public AnalysedPatientReport run(@NotNull final String runDirectory, @Nullable final String comments) throws IOException {
        final RunContext run = ProductionRunContextFactory.fromRunDirectory(runDirectory);
        assert run.isSomaticRun();

        final GenomeAnalysis genomeAnalysis = analyseGenomeData(run);
        final String tumorSample = run.tumorSample();

        final SomaticVariantAnalysis somaticVariantAnalysis = genomeAnalysis.somaticVariantAnalysis();
        final PurpleAnalysis purpleAnalysis = genomeAnalysis.purpleAnalysis();
        final StructuralVariantAnalysis structuralVariantAnalysis = genomeAnalysis.structuralVariantAnalysis();

        final List<GeneFusionData> reportableFusions = structuralVariantAnalysis.reportableFusions()
                .stream()
                .sorted(fusionComparator())
                .map(GeneFusionData::from)
                .collect(Collectors.toList());
        final List<GeneDisruptionData> reportableDisruptions = structuralVariantAnalysis.reportableDisruptions()
                .stream()
                .sorted(disruptionComparator(reporterData().panelGeneModel().transcriptMap()))
                .map(GeneDisruptionData::from)
                .collect(Collectors.toList());

        final int reportedVariantCount = somaticVariantAnalysis.variantReports().size();
        final PatientTumorLocation patientTumorLocation =
                PatientReporterHelper.extractPatientTumorLocation(baseReporterData().patientTumorLocations(), tumorSample);

        LOGGER.info("Printing analysis results:");
        LOGGER.info(" Number of variants to report : " + Integer.toString(reportedVariantCount));
        LOGGER.info("Determined copy number stats for " + Integer.toString(purpleAnalysis.genePanelSize()) + " genes which led to "
                + Integer.toString(purpleAnalysis.reportableGeneCopyNumbers().size()) + " copy numbers.");
        LOGGER.info(" Number of gene fusions to report : " + Integer.toString(reportableFusions.size()));
        LOGGER.info(" Number of gene disruptions to report : " + Integer.toString(reportableDisruptions.size()));
        LOGGER.info(" Microsatellite analysis results: " + Double.toString(somaticVariantAnalysis.indelsPerMb()) + " indels per MB");
        LOGGER.info(" Mutational load results: " + Integer.toString(somaticVariantAnalysis.mutationalLoad()));

        final Lims lims = baseReporterData().limsModel();
        final Double pathologyTumorPercentage = lims.tumorPercentageForSample(tumorSample);
        // TODO (KODU): This enrichment can be done inside variant analyser already
        final List<VariantReport> purpleEnrichedVariants = purpleAnalysis.enrichSomaticVariants(somaticVariantAnalysis.variantReports());
        final String sampleRecipient = baseReporterData().centerModel().getAddresseeStringForSample(tumorSample);

        final SampleReport sampleReport = ImmutableSampleReport.of(tumorSample,
                patientTumorLocation,
                pathologyTumorPercentage,
                lims.arrivalDateForSample(tumorSample),
                lims.arrivalDateForSample(run.refSample()),
                lims.labProceduresForSample(tumorSample),
                sampleRecipient);

        return ImmutableAnalysedPatientReport.of(sampleReport,
                purpleEnrichedVariants,
                somaticVariantAnalysis.mutationalLoad(),
                somaticVariantAnalysis.indelsPerMb(),
                purpleAnalysis.reportableGeneCopyNumbers(),
                reportableDisruptions,
                reportableFusions,
                purpleAnalysis.fittedPurity().purity(),
                purpleAnalysis.status(),
                PatientReporterHelper.findCircosPlotPath(runDirectory, tumorSample),
                Optional.ofNullable(comments),
                baseReporterData().signaturePath());
    }

    @NotNull
    private GenomeAnalysis analyseGenomeData(@NotNull RunContext run) throws IOException {
        final PurpleAnalysis purpleAnalysis = analyzePurpleCopyNumbers(run, reporterData().panelGeneModel().panel());

        final SomaticVariantAnalysis somaticVariantAnalysis = analyzeSomaticVariants(run,
                purpleAnalysis,
                variantAnalyzer(),
                reporterData().highConfidenceRegions(),
                reporterData().refGenomeFastaFile());

        final StructuralVariantAnalysis structuralVariantAnalysis = analyzeStructuralVariants(run, purpleAnalysis, structuralVariantAnalyzer());

        return ImmutableGenomeAnalysis.of(purpleAnalysis, somaticVariantAnalysis, structuralVariantAnalysis);
    }

    @NotNull
    private static PurpleAnalysis analyzePurpleCopyNumbers(@NotNull RunContext run, @NotNull Set<String> genePanel) throws IOException {
        final String runDirectory = run.runDirectory();
        final String sample = run.tumorSample();

        final PurityContext purityContext = PatientReporterHelper.loadPurity(runDirectory, sample);

        final List<PurpleCopyNumber> purpleCopyNumbers = PatientReporterHelper.loadPurpleCopyNumbers(runDirectory, sample);
        LOGGER.info(" " + purpleCopyNumbers.size() + " purple copy number regions loaded for sample " + sample);

        final List<GeneCopyNumber> panelGeneCopyNumbers = PatientReporterHelper.loadPurpleGeneCopyNumbers(runDirectory, sample)
                .stream()
                .filter(geneCopyNumber -> genePanel.contains(geneCopyNumber.gene()))
                .collect(Collectors.toList());

        return ImmutablePurpleAnalysis.builder()
                .gender(purityContext.gender())
                .status(purityContext.status())
                .fittedPurity(purityContext.bestFit())
                .fittedScorePurity(purityContext.score())
                .copyNumbers(purpleCopyNumbers)
                .panelGeneCopyNumbers(panelGeneCopyNumbers)
                .build();
    }

    @NotNull
    private static SomaticVariantAnalysis analyzeSomaticVariants(@NotNull RunContext run, @NotNull PurpleAnalysis purpleAnalysis,
            @NotNull SomaticVariantAnalyzer somaticVariantAnalyzer, @NotNull Multimap<String, GenomeRegion> highConfidenceRegions,
            @NotNull IndexedFastaSequenceFile refGenomeFastaFile) throws IOException {
        final String runDirectory = run.runDirectory();
        final String sample = run.tumorSample();

        LOGGER.info("Loading somatic variants...");
        final List<SomaticVariant> variants = PatientReporterHelper.loadPassedSomaticVariants(sample, runDirectory);
        LOGGER.info(" " + variants.size() + " PASS somatic variants loaded for sample " + sample);

        LOGGER.info("Enriching somatic variants");
        final List<EnrichedSomaticVariant> enrichedSomaticVariants =
                enrich(variants, purpleAnalysis, highConfidenceRegions, refGenomeFastaFile);

        LOGGER.info("Analyzing somatic variants....");
        return somaticVariantAnalyzer.run(enrichedSomaticVariants);
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
    private static StructuralVariantAnalysis analyzeStructuralVariants(@NotNull RunContext run, @NotNull PurpleAnalysis purpleAnalysis,
            @NotNull StructuralVariantAnalyzer structuralVariantAnalyzer) throws IOException {
        final Path structuralVariantVCF = PatientReporterHelper.findStructuralVariantVCF(run.runDirectory());
        LOGGER.info("Loading structural variants...");
        final List<StructuralVariant> structuralVariants = StructuralVariantFileLoader.fromFile(structuralVariantVCF.toString(), true);

        LOGGER.info("Enriching structural variants with purple data.");
        final List<EnrichedStructuralVariant> enrichedStructuralVariants = purpleAnalysis.enrichStructuralVariants(structuralVariants);

        LOGGER.info("Analysing structural variants...");
        return structuralVariantAnalyzer.run(enrichedStructuralVariants);
    }

    @NotNull
    private static Comparator<GeneDisruption> disruptionComparator(@NotNull final Map<String, HmfGenomeRegion> transcriptMap) {
        return Comparator.comparing(GeneDisruption::linkedAnnotation, Comparator.comparing((Transcript transcript) -> {
            final HmfGenomeRegion transcriptRegion = transcriptMap.get(transcript.transcriptId());
            final long startPosition = transcript.parent().variant().start().position();
            if (startPosition >= transcriptRegion.geneStart() && startPosition <= transcriptRegion.geneEnd()) {
                return transcript.parent().variant().start();
            } else {
                return transcript.parent().variant().end();
            }
        }));
    }

    @NotNull
    private static Comparator<GeneFusion> fusionComparator() {
        return Comparator.comparing(GeneFusion::upstreamLinkedAnnotation, transcriptComparator())
                .thenComparing(GeneFusion::downstreamLinkedAnnotation, transcriptComparator());
    }

    @NotNull
    private static Comparator<Transcript> transcriptComparator() {
        return Comparator.comparing((Transcript transcript) -> transcript.parent().variant().start())
                .thenComparing((Transcript transcript) -> transcript.parent().variant().end());
    }
}
