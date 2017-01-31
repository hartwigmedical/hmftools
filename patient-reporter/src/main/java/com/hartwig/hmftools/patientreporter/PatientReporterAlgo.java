package com.hartwig.hmftools.patientreporter;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.copynumber.CopyNumber;
import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.variant.VariantConsequence;
import com.hartwig.hmftools.common.variant.vcf.VCFFileWriter;
import com.hartwig.hmftools.common.variant.vcf.VCFSomaticFile;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberAnalysis;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberAnalyzer;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberReport;
import com.hartwig.hmftools.patientreporter.report.PDFWriter;
import com.hartwig.hmftools.patientreporter.util.ConsequenceCount;
import com.hartwig.hmftools.patientreporter.variants.VariantAnalysis;
import com.hartwig.hmftools.patientreporter.variants.VariantAnalyzer;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import net.sf.dynamicreports.report.exception.DRException;

class PatientReporterAlgo {
    private static final Logger LOGGER = LogManager.getLogger(PatientReporterAlgo.class);

    @NotNull
    private final String runDirectory;
    @NotNull
    private final CpctEcrfModel cpctEcrfModel;
    @NotNull
    private final VariantAnalyzer variantAnalyzer;
    @NotNull
    private final CopyNumberAnalyzer copyNumberAnalyzer;
    @Nullable
    private final String outputDirectory;
    @Nullable
    private final PDFWriter pdfWriter;
    private final boolean batchMode;

    PatientReporterAlgo(@NotNull final String runDirectory, @NotNull final CpctEcrfModel cpctEcrfModel,
            @NotNull final VariantAnalyzer variantAnalyzer, @NotNull final CopyNumberAnalyzer copyNumberAnalyzer,
            @Nullable final String outputDirectory, @Nullable final PDFWriter pdfWriter, final boolean batchMode) {
        this.runDirectory = runDirectory;
        this.cpctEcrfModel = cpctEcrfModel;
        this.variantAnalyzer = variantAnalyzer;
        this.copyNumberAnalyzer = copyNumberAnalyzer;
        this.outputDirectory = outputDirectory;
        this.pdfWriter = pdfWriter;
        this.batchMode = batchMode;
    }

    void run() throws IOException, HartwigException, DRException {
        if (batchMode) {
            batchRun();
        } else {
            final PatientReport report = patientRun();
            if (pdfWriter != null) {
                final String pdfReport = pdfWriter.write(report);
                LOGGER.info("  Written PDF report to " + pdfReport);

            }
        }
    }

    private void batchRun() throws IOException, HartwigException {
        // KODU: We assume "run directory" is a path with a lot of directories on which we can all run in patient mode.
        final VariantConsequence[] consequences = VariantConsequence.values();

        String header = "SAMPLE,VARIANT_COUNT,PASS_ONLY_COUNT,CONSENSUS_COUNT,MISSENSE_COUNT,CONSEQUENCE_COUNT";
        for (final VariantConsequence consequence : consequences) {
            header += ("," + consequence.name() + "_COUNT");
        }
        System.out.println(header);

        for (final Path run : Files.list(new File(runDirectory).toPath()).collect(Collectors.toList())) {
            final VCFSomaticFile variantFile = PatientReporterHelper.loadVariantFile(run.toFile().getPath());
            final VariantAnalysis analysis = variantAnalyzer.run(variantFile.variants());

            final Map<VariantConsequence, Integer> counts = ConsequenceCount.count(analysis.consensusPassedVariants());
            String consequenceList = Strings.EMPTY;
            for (final VariantConsequence consequence : consequences) {
                consequenceList += ("," + counts.get(consequence).toString());
            }

            System.out.println(
                    variantFile.sample() + "," + variantFile.variants().size() + "," + analysis.passedVariants().size()
                            + "," + analysis.consensusPassedVariants().size() + ","
                            + analysis.missenseVariants().size() + "," + analysis.findings().size() + consequenceList);
        }
    }

    @NotNull
    private PatientReport patientRun() throws IOException, HartwigException {
        LOGGER.info(" Loading genomic data...");
        final VCFSomaticFile variantFile = PatientReporterHelper.loadVariantFile(runDirectory);
        final String sample = variantFile.sample();
        LOGGER.info("  " + variantFile.variants().size() + " somatic variants loaded for " + sample);

        final List<CopyNumber> copyNumbers = PatientReporterHelper.loadCNVFile(runDirectory, sample);
        LOGGER.info("  " + copyNumbers.size() + " copy number regions loaded for sample " + sample);

        LOGGER.info(" Analyzing data...");
        final VariantAnalysis variantAnalysis = variantAnalyzer.run(variantFile.variants());
        final CopyNumberAnalysis copyNumberAnalysis = copyNumberAnalyzer.run(copyNumbers);

        final int passedCount = variantAnalysis.passedVariants().size();
        final int consensusPassedCount = variantAnalysis.consensusPassedVariants().size();
        final int mutationalLoad = variantAnalysis.missenseVariants().size();
        final int consequentialVariantCount = variantAnalysis.consequentialVariants().size();
        final int potentialMNVCount = variantAnalysis.potentialConsequentialMNVs().size();

        LOGGER.info("  Number of variants after applying pass-only filter : " + Integer.toString(passedCount));
        LOGGER.info("  Number of variants after applying consensus rule : " + Integer.toString(consensusPassedCount));
        LOGGER.info("  Number of missense variants in consensus rule (mutational load) : " + Integer.toString(
                mutationalLoad));
        LOGGER.info("  Number of consequential variants to report : " + Integer.toString(consequentialVariantCount));
        LOGGER.info("  Number of potential consequential MNVs : " + Integer.toString(potentialMNVCount));
        LOGGER.info("  Determined copy number stats for " + copyNumberAnalysis.stats().size() + " regions.");

        if (outputDirectory != null) {
            writeToFiles(outputDirectory + File.separator + sample, variantAnalysis, copyNumberAnalysis);
        }

        final String tumorType = PatientReporterHelper.extractTumorType(cpctEcrfModel, sample);
        return new PatientReport(sample, variantAnalysis.findings(), copyNumberAnalysis.findings(), mutationalLoad,
                tumorType);
    }

    private static void writeToFiles(@NotNull final String baseName, @NotNull final VariantAnalysis variantAnalysis,
            @NotNull final CopyNumberAnalysis copyNumberAnalysis) throws IOException {
        final String consensusVCF = baseName + "_consensus_variants.vcf";
        VCFFileWriter.writeSomaticVCF(consensusVCF, variantAnalysis.consensusPassedVariants());
        LOGGER.info("    Written consensus-passed variants to " + consensusVCF);

        final String consequentialVCF = baseName + "_consequential_variants.vcf";
        VCFFileWriter.writeSomaticVCF(consequentialVCF, variantAnalysis.consequentialVariants());
        LOGGER.info("    Written consequential variants to " + consequentialVCF);

        final String varReportFile = baseName + "_variant_report.csv";
        final List<VariantReport> varFindings = variantAnalysis.findings();
        Files.write(new File(varReportFile).toPath(), varToCSV(varFindings));
        LOGGER.info("    Written " + varFindings.size() + " variants to report to " + varReportFile);

        final String cnvReportFile = baseName + "_copynumber_report.csv";
        final List<CopyNumberReport> cnvFindings = copyNumberAnalysis.findings();
        Files.write(new File(cnvReportFile).toPath(), cnvToCSV(cnvFindings));
        LOGGER.info("    Written " + cnvFindings.size() + " copy-numbers to report to " + cnvReportFile);
    }

    @NotNull
    private static List<String> varToCSV(@NotNull final List<VariantReport> reports) {
        final List<String> lines = Lists.newArrayList();
        lines.add("GENE,POSITION,REF,ALT,TRANSCRIPT,CDS,AA,CONSEQUENCE,COSMIC_ID,ALLELE_READ_COUNT,TOTAL_READ_COUNT");
        lines.addAll(reports.stream().map(
                report -> report.gene() + "," + report.position() + "," + report.ref() + "," + report.alt() + ","
                        + report.transcript() + "," + report.hgvsCoding() + "," + report.hgvsProtein() + ","
                        + report.consequence() + "," + report.cosmicID() + "," + Integer.toString(
                        report.alleleReadCount()) + "," + Integer.toString(report.totalReadCount())).
                collect(Collectors.toList()));
        return lines;
    }

    @NotNull
    private static List<String> cnvToCSV(@NotNull final List<CopyNumberReport> reports) {
        final List<String> lines = Lists.newArrayList();
        lines.add("GENE,TRANSCRIPT,FINDING");
        lines.addAll(reports.stream().map(report -> report.gene() + "," + report.transcript() + "," + Integer.toString(
                report.copyNumber())).collect(Collectors.toList()));
        return lines;
    }
}
