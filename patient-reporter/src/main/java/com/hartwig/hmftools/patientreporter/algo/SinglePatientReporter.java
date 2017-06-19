package com.hartwig.hmftools.patientreporter.algo;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.hartwig.hmftools.common.copynumber.CopyNumber;
import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.lims.LimsModel;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityScore;
import com.hartwig.hmftools.common.variant.vcf.VCFFileWriter;
import com.hartwig.hmftools.common.variant.vcf.VCFSomaticFile;
import com.hartwig.hmftools.patientreporter.PatientReport;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberAnalysis;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberAnalyzer;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberReport;
import com.hartwig.hmftools.patientreporter.purple.ImmutablePurpleAnalysis;
import com.hartwig.hmftools.patientreporter.purple.PurpleAnalysis;
import com.hartwig.hmftools.patientreporter.util.FindingsToCSV;
import com.hartwig.hmftools.patientreporter.util.PatientReportFormat;
import com.hartwig.hmftools.patientreporter.variants.VariantAnalysis;
import com.hartwig.hmftools.patientreporter.variants.VariantAnalyzer;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class SinglePatientReporter {
    private static final Logger LOGGER = LogManager.getLogger(SinglePatientReporter.class);

    @NotNull
    private final CpctEcrfModel cpctEcrfModel;
    @NotNull
    private final LimsModel limsModel;
    @NotNull
    private final VariantAnalyzer variantAnalyzer;
    @NotNull
    private final CopyNumberAnalyzer copyNumberAnalyzer;
    @Nullable
    private final String tmpDirectory;
    private final boolean usePurple;

    public SinglePatientReporter(boolean usePurple, @NotNull final CpctEcrfModel cpctEcrfModel, @NotNull final LimsModel limsModel,
            @NotNull final VariantAnalyzer variantAnalyzer, @NotNull final CopyNumberAnalyzer copyNumberAnalyzer,
            @Nullable final String tmpDirectory) {
        this.cpctEcrfModel = cpctEcrfModel;
        this.limsModel = limsModel;
        this.variantAnalyzer = variantAnalyzer;
        this.copyNumberAnalyzer = copyNumberAnalyzer;
        this.tmpDirectory = tmpDirectory;
        this.usePurple = usePurple;
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
        LOGGER.info("  Number of missense variants in consensus rule (mutational load) : " + Integer.toString(
                mutationalLoad));
        LOGGER.info("  Number of consequential variants to report : " + Integer.toString(consequentialVariantCount));
        LOGGER.info("  Number of potential consequential MNVs : " + Integer.toString(potentialMNVCount));
        if (potentialMNVCount > 0) {
            LOGGER.warn(" !! Non-zero number of potentials MNV ");
        }
        LOGGER.info("  Determined copy number stats for " + Integer.toString(copyNumberAnalysis.stats().size())
                + " regions which led to " + Integer.toString(copyNumberAnalysis.findings().size()) + " findings.");

        final String tumorType = PatientReporterHelper.extractTumorType(cpctEcrfModel, sample);
        final Double tumorPercentage = limsModel.findTumorPercentageForSample(sample);
        final List<VariantReport> purpleEnrichedVariants = purpleAnalysis.enrich(variantAnalysis.findings());
        return new PatientReport(sample, purpleEnrichedVariants, copyNumberAnalysis.findings(), mutationalLoad,
                tumorType, tumorPercentage, purpleAnalysis.fittedPurity());
    }

    @NotNull
    public GenomeAnalysis analyseGenomeData(@NotNull final String runDirectory) throws IOException, HartwigException {
        LOGGER.info(" Loading somatic variants...");
        final VCFSomaticFile variantFile = PatientReporterHelper.loadVariantFile(runDirectory);
        final String sample = variantFile.sample();
        LOGGER.info("  " + variantFile.variants().size() + " somatic variants loaded for sample " + sample);

        LOGGER.info(" Loading freec somatic copy numbers...");
        final List<CopyNumber> copyNumbers = PatientReporterHelper.loadCNVFile(runDirectory, sample);
        LOGGER.info("  " + copyNumbers.size() + " freec copy number regions loaded for sample " + sample);

        LOGGER.info(" Loading purity numbers...");
        final FittedPurity purity = PatientReporterHelper.loadPurity(runDirectory, sample);
        final FittedPurityScore purityScore = PatientReporterHelper.loadPurityScore(runDirectory, sample);
        final List<PurpleCopyNumber> purpleCopyNumbers = PatientReporterHelper.loadPurpleCopyNumbers(runDirectory,
                sample);
        LOGGER.info("  " + purpleCopyNumbers.size() + " purple copy number regions loaded for sample " + sample);
        final PurpleAnalysis purpleAnalysis = ImmutablePurpleAnalysis.of(purity, purityScore, purpleCopyNumbers);
        if (Doubles.greaterThan(purpleAnalysis.purityUncertainty(), 0.03)) {
            LOGGER.warn("Purity uncertainty (" + PatientReportFormat.formatPercent(purpleAnalysis.purityUncertainty())
                    + ") range exceeds 3%. Proceed with caution.");
        }

        LOGGER.info(" Analyzing somatics....");
        final VariantAnalysis variantAnalysis = variantAnalyzer.run(variantFile.variants());

        LOGGER.info(" Analyzing {} somatic copy numbers...", usePurple ? "purple" : "freec");
        final CopyNumberAnalysis copyNumberAnalysis = copyNumberAnalyzer.run(usePurple ? purpleAnalysis.ploidyAdjustedCopyNumbers() : copyNumbers);

        if (tmpDirectory != null) {
            writeIntermediateDataToTmpFiles(tmpDirectory, variantFile, variantAnalysis, copyNumberAnalysis);
        }

        return new GenomeAnalysis(sample, variantAnalysis, copyNumberAnalysis, purpleAnalysis);
    }

    private static void writeIntermediateDataToTmpFiles(@NotNull final String basePath,
            @NotNull final VCFSomaticFile originalFile, @NotNull final VariantAnalysis variantAnalysis,
            @NotNull final CopyNumberAnalysis copyNumberAnalysis) throws IOException {
        final String baseName = basePath + File.separator + originalFile.sample();
        final String consensusVCF = baseName + "_consensus_variants.vcf";
        VCFFileWriter.writeSomaticVCF(consensusVCF, originalFile, variantAnalysis.consensusPassedVariants());
        LOGGER.info("    Written consensus-passed variants to " + consensusVCF);

        final String consequentialVCF = baseName + "_consequential_variants.vcf";
        VCFFileWriter.writeSomaticVCF(consequentialVCF, originalFile, variantAnalysis.consequentialVariants());
        LOGGER.info("    Written consequential variants to " + consequentialVCF);

        final String varReportFile = baseName + "_variant_report.csv";
        final List<VariantReport> varFindings = variantAnalysis.findings();
        Files.write(new File(varReportFile).toPath(), FindingsToCSV.varToCSV(varFindings));
        LOGGER.info("    Written " + varFindings.size() + " variants to report to " + varReportFile);

        final String cnvReportFile = baseName + "_copynumber_report.csv";
        final List<CopyNumberReport> cnvFindings = copyNumberAnalysis.findings();
        Files.write(new File(cnvReportFile).toPath(), FindingsToCSV.cnvToCSV(cnvFindings));
        LOGGER.info("    Written " + cnvFindings.size() + " copy-numbers to report to " + cnvReportFile);
    }
}
