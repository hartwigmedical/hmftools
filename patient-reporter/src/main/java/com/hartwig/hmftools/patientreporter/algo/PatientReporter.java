package com.hartwig.hmftools.patientreporter.algo;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.ecrf.doid.TumorLocationDoidMapping;
import com.hartwig.hmftools.common.ecrf.projections.PatientCancerType;
import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityScore;
import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFileLoader;
import com.hartwig.hmftools.patientreporter.BaseReporterData;
import com.hartwig.hmftools.patientreporter.HmfReporterData;
import com.hartwig.hmftools.patientreporter.ImmutableSampleReport;
import com.hartwig.hmftools.patientreporter.ImmutableSequencedPatientReport;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.SequencedPatientReport;
import com.hartwig.hmftools.patientreporter.civic.AlterationAnalyzer;
import com.hartwig.hmftools.patientreporter.copynumber.ImmutablePurpleAnalysis;
import com.hartwig.hmftools.patientreporter.copynumber.PurpleAnalysis;
import com.hartwig.hmftools.patientreporter.report.data.Alteration;
import com.hartwig.hmftools.patientreporter.report.data.GeneDisruptionData;
import com.hartwig.hmftools.patientreporter.report.data.GeneFusionData;
import com.hartwig.hmftools.patientreporter.variants.VariantAnalysis;
import com.hartwig.hmftools.patientreporter.variants.VariantAnalyzer;
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
    public abstract VariantAnalyzer variantAnalyzer();

    @NotNull
    public abstract StructuralVariantAnalyzer structuralVariantAnalyzer();

    @NotNull
    public abstract AlterationAnalyzer civicAnalyzer();

    @NotNull
    public SequencedPatientReport run(@NotNull final String runDirectory, @Nullable final String comments) throws IOException {
        final RunContext run = ProductionRunContextFactory.fromRunDirectory(runDirectory);
        final GenomeAnalysis genomeAnalysis = analyseGenomeData(run.tumorSample(), runDirectory);
        assert run.isSomaticRun() && run.tumorSample().equals(genomeAnalysis.sample());

        final String tumorSample = genomeAnalysis.sample();
        final VariantAnalysis variantAnalysis = genomeAnalysis.variantAnalysis();
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

        final int passedVariantCount = variantAnalysis.passedVariants().size();
        final int mutationalLoad = variantAnalysis.mutationalLoad();
        final int consequentialVariantCount = variantAnalysis.consequentialVariants().size();
        final int structuralVariantCount = structuralVariantAnalysis.annotations().size();
        final PatientCancerType patientCancerType =
                PatientReporterHelper.extractPatientCancerType(baseReporterData().patientsCancerTypes(), tumorSample);

        final List<Alteration> alterations;
        if (patientCancerType != null) {
            final TumorLocationDoidMapping doidMapping = TumorLocationDoidMapping.fromResource("/tumor_location_doid_mapping.csv");
            alterations = civicAnalyzer().run(variantAnalysis.findings(),
                    purpleAnalysis.reportableGeneCopyNumbers(),
                    reportableDisruptions,
                    reportableFusions,
                    reporterData().panelGeneModel(),
                    doidMapping.doidsForTumorType(patientCancerType.primaryTumorLocation()));
        } else {
            LOGGER.warn("Could not run civic analyzer as (curated) primary tumor location is not known");
            alterations = Lists.newArrayList();
        }

        LOGGER.info(" Printing analysis results:");
        LOGGER.info("  Number of passed variants : " + Integer.toString(passedVariantCount));
        LOGGER.info("  Number of missense variants (mutational load) : " + Integer.toString(mutationalLoad));
        LOGGER.info("  Number of consequential variants to report : " + Integer.toString(consequentialVariantCount));
        LOGGER.info(" Determined copy number stats for " + Integer.toString(purpleAnalysis.genePanelSize()) + " genes which led to "
                + Integer.toString(purpleAnalysis.reportableGeneCopyNumbers().size()) + " copy numbers.");
        LOGGER.info("  Number of unreported structural variants : " + Integer.toString(structuralVariantCount));
        LOGGER.info("  Number of gene fusions to report : " + Integer.toString(reportableFusions.size()));
        LOGGER.info("  Number of gene disruptions to report : " + Integer.toString(reportableDisruptions.size()));
        LOGGER.info("  Number of CIViC alterations to report : " + alterations.size());
        LOGGER.info("  Microsatellite analysis results: " + variantAnalysis.indelsPerMb() + " indels per MB");

        final Lims lims = baseReporterData().limsModel();
        final Double tumorPercentage = lims.tumorPercentageForSample(tumorSample);
        final List<VariantReport> purpleEnrichedVariants = purpleAnalysis.enrichSomaticVariants(variantAnalysis.findings());
        final String sampleRecipient = baseReporterData().centerModel().getAddresseeStringForSample(tumorSample);

        final SampleReport sampleReport = ImmutableSampleReport.of(tumorSample, patientCancerType,
                tumorPercentage,
                lims.arrivalDateForSample(tumorSample),
                lims.arrivalDateForSample(run.refSample()),
                lims.labProceduresForSample(tumorSample),
                sampleRecipient);

        return ImmutableSequencedPatientReport.of(sampleReport,
                purpleEnrichedVariants,
                mutationalLoad,
                variantAnalysis.indelsPerMb(),
                purpleAnalysis.reportableGeneCopyNumbers(),
                reportableDisruptions,
                reportableFusions,
                purpleAnalysis.purityString(),
                alterations,
                PatientReporterHelper.findCircosPlotPath(runDirectory, tumorSample),
                Optional.ofNullable(comments),
                baseReporterData().signaturePath());
    }

    @NotNull
    private GenomeAnalysis analyseGenomeData(@NotNull final String sample, @NotNull final String runDirectory) throws IOException {
        LOGGER.info(" Loading somatic snv and indels...");
        final List<SomaticVariant> variants = PatientReporterHelper.loadPassedSomaticVariants(sample, runDirectory);
        LOGGER.info("  " + variants.size() + " somatic passed snps, mnps and indels loaded for sample " + sample);
        LOGGER.info(" Analyzing somatic snp/mnp and indels....");
        final VariantAnalysis variantAnalysis = variantAnalyzer().run(variants);

        LOGGER.info(" Loading purity numbers...");
        final PurityContext context = PatientReporterHelper.loadPurity(runDirectory, sample);
        if (context.status().equals(FittedPurityStatus.NO_TUMOR)) {
            LOGGER.warn("PURPLE DID NOT DETECT A TUMOR. Proceed with utmost caution!");
        }

        final FittedPurity purity = context.bestFit();
        final FittedPurityScore purityScore = context.score();
        final List<PurpleCopyNumber> purpleCopyNumbers = PatientReporterHelper.loadPurpleCopyNumbers(runDirectory, sample);
        final List<GeneCopyNumber> panelGeneCopyNumbers = PatientReporterHelper.loadPurpleGeneCopyNumbers(runDirectory, sample)
                .stream()
                .filter(x -> reporterData().panelGeneModel().panel().contains(x.gene()))
                .collect(Collectors.toList());

        LOGGER.info("  " + purpleCopyNumbers.size() + " purple copy number regions loaded for sample " + sample);
        LOGGER.info(" Analyzing purple somatic copy numbers...");
        final PurpleAnalysis purpleAnalysis = ImmutablePurpleAnalysis.builder()
                .gender(context.gender())
                .status(context.status())
                .fittedPurity(purity)
                .fittedScorePurity(purityScore)
                .copyNumbers(purpleCopyNumbers)
                .panelGeneCopyNumbers(panelGeneCopyNumbers)
                .build();

        final Path structuralVariantVCF = PatientReporterHelper.findStructuralVariantVCF(runDirectory);
        LOGGER.info(" Loading structural variants...");
        final List<StructuralVariant> structuralVariants = StructuralVariantFileLoader.fromFile(structuralVariantVCF.toString(), true);
        LOGGER.info(" Enriching structural variants with purple data.");
        final List<EnrichedStructuralVariant> enrichedStructuralVariants = purpleAnalysis.enrichStructuralVariants(structuralVariants);
        LOGGER.info(" Analysing structural variants...");
        final StructuralVariantAnalysis structuralVariantAnalysis = structuralVariantAnalyzer().run(enrichedStructuralVariants, false);
        return ImmutableGenomeAnalysis.of(sample, variantAnalysis, purpleAnalysis, structuralVariantAnalysis);
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
