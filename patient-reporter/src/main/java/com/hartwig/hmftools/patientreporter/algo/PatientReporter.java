package com.hartwig.hmftools.patientreporter.algo;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
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
import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.BaseReportData;
import com.hartwig.hmftools.patientreporter.ImmutableAnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.ImmutableSampleReport;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.SequencedReportData;
import com.hartwig.hmftools.patientreporter.copynumber.ImmutablePurpleAnalysis;
import com.hartwig.hmftools.patientreporter.copynumber.PurpleAnalysis;
import com.hartwig.hmftools.patientreporter.variants.SomaticVariantAnalysis;
import com.hartwig.hmftools.patientreporter.variants.SomaticVariantAnalyzer;
import com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalysis;
import com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalyzer;
import com.hartwig.hmftools.svannotation.annotations.GeneDisruption;
import com.hartwig.hmftools.svannotation.annotations.GeneFusion;

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
    public abstract BaseReportData baseReportData();

    @NotNull
    public abstract SequencedReportData sequencedReportData();

    @NotNull
    public abstract StructuralVariantAnalyzer structuralVariantAnalyzer();

    @NotNull
    public AnalysedPatientReport run(@NotNull final String runDirectory, @Nullable final String comments) throws IOException {
        final RunContext run = ProductionRunContextFactory.fromRunDirectory(runDirectory);
        assert run.isSomaticRun();

        final PurpleAnalysis purpleAnalysis = analyzePurpleCopyNumbers(run, sequencedReportData().panelGeneModel());
        final List<GeneCopyNumber> reportableGeneCopynumbers =
                purpleAnalysis.reportableGeneCopyNumbers(sequencedReportData().panelGeneModel());

        final String tumorSample = run.tumorSample();
        PatientTumorLocation patientTumorLocation =
                PatientReporterFileLoader.extractPatientTumorLocation(baseReportData().patientTumorLocations(), tumorSample);

        final StructuralVariantAnalysis structuralVariantAnalysis =
                analyzeStructuralVariants(run, purpleAnalysis, structuralVariantAnalyzer());

        final SomaticVariantAnalysis somaticVariantAnalysis = analyzeSomaticVariants(run,
                purpleAnalysis,
                sequencedReportData().somaticVariantEnrichment(),
                sequencedReportData().panelGeneModel(),
                sequencedReportData().highConfidenceRegions(),
                sequencedReportData().refGenomeFastaFile(),
                patientTumorLocation,
                structuralVariantAnalysis.fusions(),
                sequencedReportData().actionabilityAnalyzer());

        final List<GeneFusion> reportableFusions = structuralVariantAnalysis.reportableFusions();
        final List<GeneDisruption> reportableDisruptions = structuralVariantAnalysis.reportableDisruptions();

        LOGGER.info("Printing analysis results:");
        LOGGER.info(" Number of somatic variants to report : " + Integer.toString(somaticVariantAnalysis.variantsToReport().size()));
        LOGGER.info(" Microsatellite analysis results: " + Double.toString(somaticVariantAnalysis.microsatelliteIndelsPerMb())
                + " indels per MB");
        LOGGER.info(" Mutational load results: " + Integer.toString(somaticVariantAnalysis.tumorMutationalLoad()));
        LOGGER.info(" Tumor mutational burden: " + Double.toString(somaticVariantAnalysis.tumorMutationalBurden())
                + " number of mutations per MB");
        LOGGER.info(" Number of copy number events to report: " + Integer.toString(reportableGeneCopynumbers.size()));
        LOGGER.info(" Number of gene fusions to report : " + Integer.toString(reportableFusions.size()));
        LOGGER.info(" Number of gene disruptions to report : " + Integer.toString(reportableDisruptions.size()));

        LOGGER.info("Printing evidence results:");
        LOGGER.info(" Number of actionability variants to report: " + Integer.toString(somaticVariantAnalysis.evidencePerVariant().size()));
        LOGGER.info(" Number of actionability gene copy numbers to report: " + Integer.toString(somaticVariantAnalysis.evidencePerCopyNumber().size()));
       // LOGGER.info("Number of actionability fusions to report: " + Integer.toString(somaticVariantAnalysis.evidenceFusions().size()));

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

        for (Map.Entry<GeneCopyNumber, List<EvidenceItem>> evidencePerCNV : somaticVariantAnalysis.evidencePerCopyNumber().entrySet()) {
            evidenceItems.addAll(evidencePerCNV.getValue());
        }

        return ImmutableAnalysedPatientReport.of(sampleReport,
                purpleAnalysis.status(),
                purpleAnalysis.fittedPurity().purity(),
                evidenceItems,
                somaticVariantAnalysis.variantsToReport(),
                somaticVariantAnalysis.driverCatalog(),
                somaticVariantAnalysis.microsatelliteIndelsPerMb(),
                somaticVariantAnalysis.tumorMutationalLoad(),
                somaticVariantAnalysis.tumorMutationalBurden(),
                Lists.newArrayList(),
                reportableGeneCopynumbers,
                reportableFusions,
                reportableDisruptions,
                PatientReporterFileLoader.findCircosPlotPath(runDirectory, tumorSample),
                Optional.ofNullable(comments),
                baseReportData().signaturePath());
    }

    @NotNull
    private static PurpleAnalysis analyzePurpleCopyNumbers(@NotNull RunContext run, @NotNull GeneModel panelGeneModel) throws IOException {
        final String runDirectory = run.runDirectory();
        final String sample = run.tumorSample();

        LOGGER.info("Loading purple data for sample " + sample);
        final PurityContext purityContext = PatientReporterFileLoader.loadPurity(runDirectory, sample);

        final List<PurpleCopyNumber> purpleCopyNumbers = PatientReporterFileLoader.loadPurpleCopyNumbers(runDirectory, sample);
        LOGGER.info(" " + purpleCopyNumbers.size() + " purple copy number regions loaded for sample " + sample);

        final List<GeneCopyNumber> panelGeneCopyNumbers = PatientReporterFileLoader.loadPurpleGeneCopyNumbers(runDirectory, sample)
                .stream()
                .filter(geneCopyNumber -> panelGeneModel.cnvGenePanel().contains(geneCopyNumber.gene()))
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
            @NotNull SomaticEnrichment somaticEnrichment, @NotNull GeneModel geneModel,
            @NotNull Multimap<String, GenomeRegion> highConfidenceRegions, @NotNull IndexedFastaSequenceFile refGenomeFastaFile,
            @Nullable PatientTumorLocation patientTumorLocation, @NotNull List<GeneFusion> fusions,
            @NotNull ActionabilityAnalyzer actionabilityAnalyzerData) throws IOException {
        final String runDirectory = run.runDirectory();
        final String sample = run.tumorSample();

        LOGGER.info("Loading somatic variants...");
        final List<SomaticVariant> variants = PatientReporterFileLoader.loadPassedSomaticVariants(sample, runDirectory, somaticEnrichment);
        LOGGER.info(" " + variants.size() + " PASS somatic variants loaded for sample " + sample);

        LOGGER.info("Enriching somatic variants");
        final List<EnrichedSomaticVariant> enrichedSomaticVariants =
                enrich(variants, purpleAnalysis, highConfidenceRegions, refGenomeFastaFile);

        LOGGER.info("Analyzing somatic variants....");
        return SomaticVariantAnalyzer.run(enrichedSomaticVariants,
                geneModel.somaticVariantGenePanel(),
                geneModel.geneDriverCategoryMap(),
                patientTumorLocation,
                PatientReporterFileLoader.loadPurpleGeneCopyNumbers(runDirectory, sample),
                fusions,
                actionabilityAnalyzerData, purpleAnalysis.fittedPurity().ploidy());
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
        final Path structuralVariantVCF = PatientReporterFileLoader.findStructuralVariantVCF(run.runDirectory());
        LOGGER.info("Loading structural variants...");
        final List<StructuralVariant> structuralVariants = StructuralVariantFileLoader.fromFile(structuralVariantVCF.toString(), true);

        LOGGER.info("Enriching structural variants with purple data.");
        final PurityAdjuster purityAdjuster = new PurityAdjuster(purpleAnalysis.gender(), purpleAnalysis.fittedPurity());
        final Multimap<String, PurpleCopyNumber> copyNumberMap =
                Multimaps.index(purpleAnalysis.copyNumbers(), PurpleCopyNumber::chromosome);

        final List<EnrichedStructuralVariant> enrichedStructuralVariants =
                EnrichedStructuralVariantFactory.enrich(structuralVariants, purityAdjuster, copyNumberMap);

        LOGGER.info("Analysing structural variants...");
        return structuralVariantAnalyzer.run(enrichedStructuralVariants);
    }
}
