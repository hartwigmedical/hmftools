package com.hartwig.hmftools.patientreporter.algo;

import java.io.IOException;
import java.nio.file.Path;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.ecrf.doid.TumorLocationDoidMapping;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityScore;
import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFileLoader;
import com.hartwig.hmftools.patientreporter.BaseReporterData;
import com.hartwig.hmftools.patientreporter.HmfReporterData;
import com.hartwig.hmftools.patientreporter.ImmutableSampleReport;
import com.hartwig.hmftools.patientreporter.ImmutableSequencedPatientReport;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.SequencedPatientReport;
import com.hartwig.hmftools.patientreporter.civic.CivicAnalysis;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberAnalysis;
import com.hartwig.hmftools.patientreporter.purple.ImmutablePurpleAnalysis;
import com.hartwig.hmftools.patientreporter.purple.PurpleAnalysis;
import com.hartwig.hmftools.patientreporter.report.data.Alteration;
import com.hartwig.hmftools.patientreporter.report.data.GeneDisruptionData;
import com.hartwig.hmftools.patientreporter.report.data.GeneFusionData;
import com.hartwig.hmftools.patientreporter.util.PatientReportFormat;
import com.hartwig.hmftools.patientreporter.variants.VariantAnalysis;
import com.hartwig.hmftools.patientreporter.variants.VariantAnalyzer;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;
import com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalysis;
import com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalyzer;

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
    public SequencedPatientReport run(@NotNull final String runDirectory, @Nullable final String comments)
            throws IOException, HartwigException {
        final RunContext run = ProductionRunContextFactory.fromRunDirectory(runDirectory);
        final GenomeAnalysis genomeAnalysis = analyseGenomeData(run.tumorSample(), runDirectory);
        assert run.isSomaticRun() && run.tumorSample().equals(genomeAnalysis.sample());

        final String tumorSample = genomeAnalysis.sample();
        final VariantAnalysis variantAnalysis = genomeAnalysis.variantAnalysis();
        final CopyNumberAnalysis copyNumberAnalysis = genomeAnalysis.purpleAnalysis().copyNumberAnalysis();
        final PurpleAnalysis purpleAnalysis = genomeAnalysis.purpleAnalysis();
        final StructuralVariantAnalysis svAnalysis = genomeAnalysis.structuralVariantAnalysis();
        final List<GeneFusionData> reportableFusions =
                svAnalysis.reportableFusions().stream().map(GeneFusionData::from).collect(Collectors.toList());
        final List<GeneDisruptionData> reportableDisruptions =
                svAnalysis.reportableDisruptions().stream().map(GeneDisruptionData::from).collect(Collectors.toList());

        final int totalVariantCount = variantAnalysis.allVariants().size();
        final int passedVariantCount = variantAnalysis.passedVariants().size();
        final int mutationalLoad = variantAnalysis.mutationalLoad();
        final int consequentialVariantCount = variantAnalysis.consequentialVariants().size();
        final int potentialMNVCount = variantAnalysis.potentialConsequentialMNVs().size();
        final int svCount = svAnalysis.annotations().size();
        final String tumorType = PatientReporterHelper.extractTumorType(baseReporterData().cpctEcrfModel(), tumorSample);

        final TumorLocationDoidMapping doidMapping = TumorLocationDoidMapping.fromResource("/tumor_location_doid_mapping.csv");
        final List<Alteration> alterations = CivicAnalysis.run(variantAnalysis.findings(),
                copyNumberAnalysis.findings(),
                reporterData().geneModel(),
                doidMapping.doidsForTumorType(tumorType));

        LOGGER.info(" Printing analysis results:");
        LOGGER.info("  Number of variants: " + Integer.toString(totalVariantCount));
        LOGGER.info("  Number of variants after applying pass-only filter : " + Integer.toString(passedVariantCount));
        LOGGER.info("  Number of missense variants (mutational load) : " + Integer.toString(mutationalLoad));
        LOGGER.info("  Number of consequential variants to report : " + Integer.toString(consequentialVariantCount));
        LOGGER.info("  Number of potential consequential MNVs : " + Integer.toString(potentialMNVCount));
        if (potentialMNVCount > 0) {
            LOGGER.warn(" !! Non-zero number of potentials MNV ");
            LOGGER.warn(variantAnalysis.potentialConsequentialMNVs());
        }
        LOGGER.info("  Determined copy number stats for " + Integer.toString(copyNumberAnalysis.genePanelSize()) + " genes which led to "
                + Integer.toString(copyNumberAnalysis.findings().size()) + " findings.");
        LOGGER.info("  Number of unreported structural variants : " + Integer.toString(svCount));
        LOGGER.info("  Number of gene fusions to report : " + Integer.toString(reportableFusions.size()));
        LOGGER.info("  Number of gene disruptions to report : " + Integer.toString(reportableDisruptions.size()));
        LOGGER.info("  Number of CIViC alterations to report : " + alterations.size());

        final Lims lims = baseReporterData().limsModel();
        final Double tumorPercentage = lims.tumorPercentageForSample(tumorSample);
        final List<VariantReport> purpleEnrichedVariants = purpleAnalysis.enrich(variantAnalysis.findings());
        final String sampleRecipient = baseReporterData().centerModel().getAddresseeStringForSample(tumorSample);
        final SampleReport sampleReport = ImmutableSampleReport.of(tumorSample,
                tumorType,
                tumorPercentage,
                lims.arrivalDateForSample(tumorSample),
                lims.arrivalDateForSample(run.refSample()),
                lims.labProceduresForSample(tumorSample),
                sampleRecipient);
        return ImmutableSequencedPatientReport.of(sampleReport,
                purpleEnrichedVariants,
                reportableFusions,
                reportableDisruptions,
                copyNumberAnalysis.findings(),
                mutationalLoad,
                purpleAnalysis.purityString(),
                alterations,
                Optional.ofNullable(comments),
                baseReporterData().signaturePath());
    }

    @NotNull
    private GenomeAnalysis analyseGenomeData(@NotNull final String sample, @NotNull final String runDirectory)
            throws IOException, HartwigException {
        LOGGER.info(" Loading somatic snv and indels...");
        final List<SomaticVariant> variants = PatientReporterHelper.loadSomaticSNVFile(sample, runDirectory);
        LOGGER.info("  " + variants.size() + " somatic snv and indels loaded for sample " + sample);
        LOGGER.info(" Analyzing somatic snv and indels....");
        final VariantAnalysis variantAnalysis = variantAnalyzer().run(variants);

        LOGGER.info(" Loading purity numbers...");
        final PurityContext context = PatientReporterHelper.loadPurity(runDirectory, sample);
        if (context.status().equals(FittedPurityStatus.NO_TUMOR)) {
            LOGGER.warn("PURPLE DID NOT DETECT A TUMOR. Proceed with utmost caution!");
        }

        final FittedPurity purity = context.bestFit();
        final FittedPurityScore purityScore = context.score();
        final List<PurpleCopyNumber> purpleCopyNumbers = PatientReporterHelper.loadPurpleCopyNumbers(runDirectory, sample);
        final List<GeneCopyNumber> geneCopyNumbers = PatientReporterHelper.loadPurpleGeneCopyNumbers(runDirectory, sample)
                .stream()
                .filter(x -> reporterData().geneModel().panel().contains(x.gene()))
                .collect(Collectors.toList());

        LOGGER.info("  " + purpleCopyNumbers.size() + " purple copy number regions loaded for sample " + sample);
        LOGGER.info(" Analyzing purple somatic copy numbers...");
        final PurpleAnalysis purpleAnalysis = ImmutablePurpleAnalysis.builder()
                .gender(context.gender())
                .status(context.status())
                .fittedPurity(purity)
                .fittedScorePurity(purityScore)
                .copyNumbers(purpleCopyNumbers)
                .geneCopyNumbers(geneCopyNumbers)
                .build();
        if (Doubles.greaterThan(purpleAnalysis.purityUncertainty(), 0.05)) {
            LOGGER.warn("Purity uncertainty (" + PatientReportFormat.formatPercent(purpleAnalysis.purityUncertainty())
                    + ") range exceeds 5%. Proceed with caution.");
        }

        final Path svVcfPath = PatientReporterHelper.findStructuralVariantVCF(runDirectory);
        LOGGER.info(" Loading structural variants...");
        final List<StructuralVariant> structuralVariants = StructuralVariantFileLoader.fromFile(svVcfPath.toString());
        LOGGER.info(" Analysing structural variants...");
        final StructuralVariantAnalysis svAnalysis = structuralVariantAnalyzer().run(structuralVariants);

        return new GenomeAnalysis(sample, variantAnalysis, purpleAnalysis, svAnalysis);
    }
}
