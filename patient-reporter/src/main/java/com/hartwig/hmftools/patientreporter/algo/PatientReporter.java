package com.hartwig.hmftools.patientreporter.algo;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.lims.LimsJsonModel;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityScore;
import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.variant.vcf.VCFSomaticFile;
import com.hartwig.hmftools.patientreporter.ImmutablePatientReport;
import com.hartwig.hmftools.patientreporter.PatientReport;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberAnalysis;
import com.hartwig.hmftools.patientreporter.purple.ImmutablePurpleAnalysis;
import com.hartwig.hmftools.patientreporter.purple.PurpleAnalysis;
import com.hartwig.hmftools.patientreporter.util.PatientReportFormat;
import com.hartwig.hmftools.patientreporter.variants.VariantAnalysis;
import com.hartwig.hmftools.patientreporter.variants.VariantAnalyzer;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class PatientReporter {
    private static final Logger LOGGER = LogManager.getLogger(PatientReporter.class);

    @NotNull
    private final CpctEcrfModel cpctEcrfModel;
    @NotNull
    private final LimsJsonModel limsModel;
    @NotNull
    private final VariantAnalyzer variantAnalyzer;

    public PatientReporter(@NotNull final CpctEcrfModel cpctEcrfModel, @NotNull final LimsJsonModel limsModel,
            @NotNull final VariantAnalyzer variantAnalyzer) {
        this.cpctEcrfModel = cpctEcrfModel;
        this.limsModel = limsModel;
        this.variantAnalyzer = variantAnalyzer;
    }

    @NotNull
    public PatientReport run(@NotNull final String runDirectory) throws IOException, HartwigException {
        final GenomeAnalysis genomeAnalysis = analyseGenomeData(runDirectory);

        final String sample = genomeAnalysis.sample();
        final VariantAnalysis variantAnalysis = genomeAnalysis.variantAnalysis();
        final CopyNumberAnalysis copyNumberAnalysis = genomeAnalysis.copyNumberAnalysis();
        final PurpleAnalysis purpleAnalysis = genomeAnalysis.purpleAnalysis();

        final int passedCount = variantAnalysis.passedVariants().size();
        final int consensusPassedCount = variantAnalysis.consensusPassedVariants().size();
        final int mutationalLoad = variantAnalysis.mutationalLoad();
        final int consequentialVariantCount = variantAnalysis.consequentialVariants().size();
        final int potentialMNVCount = variantAnalysis.potentialConsequentialMNVs().size();

        LOGGER.info(" Printing analysis results:");
        LOGGER.info("  Number of variants after applying pass-only filter : " + Integer.toString(passedCount));
        LOGGER.info("  Number of variants after applying consensus rule : " + Integer.toString(consensusPassedCount));
        LOGGER.info("  Number of missense variants in consensus rule (mutational load) : " + Integer.toString(mutationalLoad));
        LOGGER.info("  Number of consequential variants to report : " + Integer.toString(consequentialVariantCount));
        LOGGER.info("  Number of potential consequential MNVs : " + Integer.toString(potentialMNVCount));
        if (potentialMNVCount > 0) {
            LOGGER.warn(" !! Non-zero number of potentials MNV ");
            LOGGER.warn(variantAnalysis.potentialConsequentialMNVs());
        }
        LOGGER.info("  Determined copy number stats for " + Integer.toString(copyNumberAnalysis.genePanelSize()) + " genes which led to "
                + Integer.toString(copyNumberAnalysis.findings().size()) + " findings.");

        final String tumorType = PatientReporterHelper.extractTumorType(cpctEcrfModel, sample);
        final Double tumorPercentage = limsModel.tumorPercentageForSample(sample);
        final List<VariantReport> purpleEnrichedVariants = purpleAnalysis.enrich(variantAnalysis.findings());
        return ImmutablePatientReport.of(sample, purpleEnrichedVariants, copyNumberAnalysis.findings(), mutationalLoad, tumorType,
                tumorPercentage, purpleAnalysis.purityString());
    }

    @NotNull
    private GenomeAnalysis analyseGenomeData(@NotNull final String runDirectory) throws IOException, HartwigException {
        LOGGER.info(" Loading somatic variants...");
        final VCFSomaticFile variantFile = PatientReporterHelper.loadVariantFile(runDirectory);
        final String sample = variantFile.sample();
        LOGGER.info("  " + variantFile.variants().size() + " somatic variants loaded for sample " + sample);

        LOGGER.info(" Loading purity numbers...");
        final PurityContext context = PatientReporterHelper.loadPurity(runDirectory, sample);
        if (context.status().equals(FittedPurityStatus.NO_TUMOR)) {
            LOGGER.warn("PURPLE DID NOT DETECT A TUMOR. Proceed with utmost caution!");
        }

        final FittedPurity purity = context.bestFit();
        final FittedPurityScore purityScore = context.score();
        final List<PurpleCopyNumber> purpleCopyNumbers = PatientReporterHelper.loadPurpleCopyNumbers(runDirectory, sample);
        final List<GeneCopyNumber> geneCopyNumbers = PatientReporterHelper.loadPurpleGeneCopyNumbers(runDirectory, sample);
        LOGGER.info("  " + purpleCopyNumbers.size() + " purple copy number regions loaded for sample " + sample);
        LOGGER.info(" Analyzing purple somatic copy numbers...");
        final PurpleAnalysis purpleAnalysis =
                ImmutablePurpleAnalysis.of(context.status(), purity, purityScore, purpleCopyNumbers, geneCopyNumbers);
        if (Doubles.greaterThan(purpleAnalysis.purityUncertainty(), 0.05)) {
            LOGGER.warn("Purity uncertainty (" + PatientReportFormat.formatPercent(purpleAnalysis.purityUncertainty())
                    + ") range exceeds 5%. Proceed with caution.");
        }

        final CopyNumberAnalysis copyNumberAnalysis = purpleAnalysis.copyNumberAnalysis();

        LOGGER.info(" Analyzing somatics....");
        final VariantAnalysis variantAnalysis = variantAnalyzer.run(variantFile.variants());

        return new GenomeAnalysis(sample, variantAnalysis, copyNumberAnalysis, purpleAnalysis);
    }
}
