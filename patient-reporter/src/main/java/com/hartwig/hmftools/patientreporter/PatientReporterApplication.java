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
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.variant.VariantConsequence;
import com.hartwig.hmftools.common.variant.vcf.VCFFileWriter;
import com.hartwig.hmftools.common.variant.vcf.VCFSomaticFile;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberAnalysis;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberAnalyzer;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberReport;
import com.hartwig.hmftools.patientreporter.report.PDFWriter;
import com.hartwig.hmftools.patientreporter.slicing.Slicer;
import com.hartwig.hmftools.patientreporter.slicing.SlicerFactory;
import com.hartwig.hmftools.patientreporter.util.ConsequenceCount;
import com.hartwig.hmftools.patientreporter.variants.VariantAnalysis;
import com.hartwig.hmftools.patientreporter.variants.VariantAnalyzer;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class PatientReporterApplication {

    private static final Logger LOGGER = LogManager.getLogger(PatientReporterApplication.class);

    private static final String RUN_DIRECTORY_ARGS_DESC = "A path towards a single rundir.";
    private static final String RUN_DIRECTORY = "rundir";

    private static final String CPCT_SLICING_BED_ARGS_DESC = "A path towards the CPCT slicing bed.";
    private static final String CPCT_SLICING_BED = "cpct_slicing_bed";

    private static final String HIGH_CONFIDENCE_BED_ARGS_DESC = "A path towards the high confidence bed.";
    private static final String HIGH_CONFIDENCE_BED = "high_confidence_bed";

    private static final String HMF_SLICING_BED_ARGS_DESC = "A path towards the HMF slicing bed.";
    private static final String HMF_SLICING_BED = "hmf_slicing_bed";

    private static final String HMF_LOGO_ARGS_DESC = "A path to the HMF logo, used in the PDF report.";
    private static final String HMF_LOGO = "hmf_logo";

    private static final String OUTPUT_DIR_ARGS_DESC = "A path where, if provided, output files will be written to.";
    private static final String OUTPUT_DIR = "output_dir";

    private static final String BATCH_MODE_ARGS_DESC = "If set, runs in batch mode (Caution!!! Korneel Only)";
    private static final String BATCH_MODE = "batch_mode";

    @NotNull
    private final String runDirectory;
    @NotNull
    private final VariantAnalyzer variantAnalyzer;
    @NotNull
    private final CopyNumberAnalyzer copyNumberAnalyzer;
    @Nullable
    private final String outputDirectory;
    @Nullable
    private final PDFWriter pdfWriter;
    private final boolean batchMode;

    public static void main(final String... args) throws ParseException, IOException, HartwigException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        final String runDir = cmd.getOptionValue(RUN_DIRECTORY);
        final String cpctSlicingBed = cmd.getOptionValue(CPCT_SLICING_BED);
        final String highConfidenceBed = cmd.getOptionValue(HIGH_CONFIDENCE_BED);
        final String hmfSlicingBed = cmd.getOptionValue(HMF_SLICING_BED);
        final String outputDirectory = cmd.getOptionValue(OUTPUT_DIR);
        final String hmfLogo = cmd.getOptionValue(HMF_LOGO);
        final boolean batchMode = cmd.hasOption(BATCH_MODE);

        if (runDir == null || cpctSlicingBed == null || highConfidenceBed == null || hmfSlicingBed == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Patient-Reporter", options);
            System.exit(1);
        }

        if (outputDirectory != null) {
            final Path outputPath = new File(outputDirectory).toPath();
            if (!Files.exists(outputPath) || !Files.isDirectory(outputPath)) {
                LOGGER.warn(OUTPUT_DIR + " has to be an existing directory!");
                System.exit(1);
            }
        }

        if (hmfLogo != null) {
            final Path logoPath = new File(hmfLogo).toPath();
            if (!Files.exists(logoPath)) {
                LOGGER.warn(HMF_LOGO + " has to be an existing file!");
                System.exit(1);
            }
        }

        PDFWriter pdfWriter = null;
        if (outputDirectory != null && hmfLogo != null) {
            pdfWriter = new PDFWriter(outputDirectory, hmfLogo);
        }

        final Slicer hmfSlicingRegion = SlicerFactory.fromBedFile(hmfSlicingBed);
        final VariantAnalyzer variantAnalyzer = VariantAnalyzer.fromSlicingRegions(hmfSlicingRegion,
                SlicerFactory.fromBedFile(highConfidenceBed), SlicerFactory.fromBedFile(cpctSlicingBed));
        final CopyNumberAnalyzer copyNumberAnalyzer = CopyNumberAnalyzer.fromHmfSlicingRegion(hmfSlicingRegion);
        new PatientReporterApplication(runDir, variantAnalyzer, copyNumberAnalyzer, outputDirectory, pdfWriter,
                batchMode).run();
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();

        options.addOption(RUN_DIRECTORY, true, RUN_DIRECTORY_ARGS_DESC);
        options.addOption(CPCT_SLICING_BED, true, CPCT_SLICING_BED_ARGS_DESC);
        options.addOption(HIGH_CONFIDENCE_BED, true, HIGH_CONFIDENCE_BED_ARGS_DESC);
        options.addOption(HMF_SLICING_BED, true, HMF_SLICING_BED_ARGS_DESC);
        options.addOption(HMF_LOGO, true, HMF_LOGO_ARGS_DESC);
        options.addOption(OUTPUT_DIR, true, OUTPUT_DIR_ARGS_DESC);
        options.addOption(BATCH_MODE, false, BATCH_MODE_ARGS_DESC);

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args)
            throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    PatientReporterApplication(@NotNull final String runDirectory, @NotNull final VariantAnalyzer variantAnalyzer,
            @NotNull final CopyNumberAnalyzer copyNumberAnalyzer, @Nullable final String outputDirectory,
            @Nullable final PDFWriter pdfWriter, final boolean batchMode) {
        this.runDirectory = runDirectory;
        this.variantAnalyzer = variantAnalyzer;
        this.copyNumberAnalyzer = copyNumberAnalyzer;
        this.outputDirectory = outputDirectory;
        this.pdfWriter = pdfWriter;
        this.batchMode = batchMode;
    }

    void run() throws IOException, HartwigException {
        if (batchMode) {
            batchRun();
        } else {
            final PatientReport report = patientRun();
            if (pdfWriter != null) {
                pdfWriter.write(report);
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
        LOGGER.info("Running patient reporter on " + runDirectory);

        LOGGER.info(" Loading data...");
        final VCFSomaticFile variantFile = PatientReporterHelper.loadVariantFile(runDirectory);
        LOGGER.info("  Somatic variants loaded : " + variantFile.variants().size());

        final List<CopyNumber> copyNumbers = PatientReporterHelper.loadCNVFile(runDirectory, variantFile.sample());
        LOGGER.info("  CNV data loaded for sample " + variantFile.sample());

        LOGGER.info(" Analyzing data...");
        final VariantAnalysis variantAnalysis = variantAnalyzer.run(variantFile.variants());
        final CopyNumberAnalysis copyNumberAnalysis = copyNumberAnalyzer.run(copyNumbers);

        LOGGER.info(
                "  Number of variants after applying pass-only filter : " + variantAnalysis.passedVariants().size());
        LOGGER.info("  Number of variants after applying consensus rule : "
                + variantAnalysis.consensusPassedVariants().size());
        LOGGER.info("  Number of missense variants in consensus rule (mutational load) : "
                + variantAnalysis.missenseVariants().size());
        LOGGER.info("  Number of consequential variants to report : " + variantAnalysis.findings().size());
        LOGGER.info("  Determined copy number stats for " + copyNumberAnalysis.stats().size() + " regions.");

        if (outputDirectory != null) {
            writeToFiles(outputDirectory + File.separator + variantFile.sample(), variantAnalysis, copyNumberAnalysis);
        }

        return new PatientReport(variantFile.sample(), variantAnalysis.findings(), copyNumberAnalysis.findings(),
                variantAnalysis.missenseVariants().size());
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
